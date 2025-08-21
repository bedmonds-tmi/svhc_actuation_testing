#include <zephyr/kernel.h>
#include <zephyr/timing/timing.h>
#include <tmi/api/pressure.h>

#include <zephyr/sys/ring_buffer.h>

#include "audio_settings.h"

static const struct device * pdm_dev = DEVICE_DT_GET(DT_NODELABEL(dmic_dev));
static const struct device * pressure_dev = DEVICE_DT_GET(DT_NODELABEL(pressure));

// -----------------------------------
// PDM microphone
// -----------------------------------
K_MEM_SLAB_DEFINE_STATIC(mem_slab, BLOCK_SIZE, BLOCK_COUNT, 4);
RING_BUF_ITEM_DECLARE_POW2(pdm_ring_buf, 15);

static int pdm_init(const struct device* dev) {
	int ret;

	struct pcm_stream_cfg stream = {
		.pcm_rate = SAMPLE_FREQUENCY,
		.pcm_width = SAMPLE_BIT_WIDTH,
		.block_size = BLOCK_SIZE,
		.mem_slab  = &mem_slab
	};

	struct dmic_cfg cfg = {
		.io = {
			/* These fields can be used to limit the PDM clock
				* configurations that the driver is allowed to use
				* to those supported by the microphone.
				*/
			.min_pdm_clk_freq = 1000000,
			.max_pdm_clk_freq = 4800000,
			.min_pdm_clk_dc   = 40, // Clk duty cycle, min: 40%, max: 60%
			.max_pdm_clk_dc   = 60
		},
		.streams = &stream,
		.channel = {
			.req_num_chan = 2,
			.req_num_streams = 1
		},
	};

	cfg.channel.req_chan_map_lo = dmic_build_channel_map(1, 0, PDM_CHAN_LEFT) | dmic_build_channel_map(0, 0, PDM_CHAN_RIGHT);

	printk("PCM output rate: %u, channels: %u", cfg.streams[0].pcm_rate, cfg.channel.req_num_chan);

	ret = dmic_configure(dev, &cfg);
	if (ret < 0) {
		printk("Failed to configure the driver: %d", ret);
		return ret;
	}

	ret = dmic_trigger(dev, DMIC_TRIGGER_START);
	if (ret < 0) {
		printk("START trigger failed: %d", ret);
		return ret;
	}

	return 0;
}

static int pdm_stop(const struct device* dev){
	int ret;
	ret = dmic_trigger(dev, DMIC_TRIGGER_STOP);
	return ret;
}

// -----------------------------------
// CF Sensor
// -----------------------------------

#define MY_STACK_SIZE 32000
#define MY_PRIORITY -2

K_THREAD_STACK_DEFINE(my_stack_area, MY_STACK_SIZE);

struct k_work_q my_work_q;

// Data management
float pressure_buf[16];
int count = 0;

extern int pressure_window_length;
extern int audio_window_length;
extern void predict_actuation(const int16_t* audio_data, const int audio_window_size, const float* pressure_data);

void buffer_pressure() {
	float pressure;
	//int64_t us = k_uptime_ticks() * 1000000LL / CONFIG_SYS_CLOCK_TICKS_PER_SEC;
	pressure_get(pressure_dev, &pressure);
	pressure_buf[count] = pressure;
	count++;
	// printk("%f\n", (double)pressure);

	if (count == 16){
		printk("Enough pressure data collected\n");
		// Make a copy of pressure buffer
		float tmp_pressure_buf[pressure_window_length];
		for(int i=0; i < 15; ++i){
			tmp_pressure_buf[i] = pressure_buf[i];
		}
		// Reset buf counter
		count = 0;

		printk("Retrieve audio data\n");
		// Get audio data window
		int16_t aud_buf[audio_window_length];		// NOTE: this is slightly confusing, but audio_window_length is the number of samples per window AFTER down-sampling
		int aud_bytes = ring_buf_get(&pdm_ring_buf, (uint8_t*)aud_buf, sizeof(aud_buf));

		printk("Start prediction function\n");
		// Create a thread to handle actuation prediction
		predict_actuation(aud_buf, (int)(aud_bytes/sizeof(int16_t)), tmp_pressure_buf);
	}
}

K_WORK_DEFINE(pressure_work, buffer_pressure);

void pressure_timer_handler(struct k_timer *dummy) {
    k_work_submit_to_queue(&my_work_q, &pressure_work);
}

K_TIMER_DEFINE(pressure_timer, pressure_timer_handler, NULL);

// -----------------------------------
// PDM mic audio collection thread management
// -----------------------------------
k_tid_t audio_tid = NULL;
struct k_thread audio_thread_data;

#define AUDIO_STACK_SIZE 5000
#define AUDIO_PRIORITY   -2

K_THREAD_STACK_DEFINE(audio_stack_area, AUDIO_STACK_SIZE);

/**
 * @brief Audio collection storage loop
 */
static void ac_subsys(const struct device* dev){
	int ret;

	while(1) {
		/* Read data from PDM microphone */
		uint8_t *buffer = NULL;
		uint32_t size = 0;

		ret = dmic_read(dev, 0, (void *)&buffer, &size, READ_TIMEOUT);
		if (ret < 0) {
			printk("read failed: %d\n", ret);
			return;
		}

		/* Copy data to ring buffer */
		ret = ring_buf_put(&pdm_ring_buf, buffer, size); // Returns number of bytes copied to ring buffer
		k_mem_slab_free(&mem_slab, (void *)buffer);
	}

	ret = pdm_stop(dev);
	if(ret < 0) {
		return;
	}
}

int ac_init_thread(const struct device* dev){
	int ret;
	ret = pdm_init(dev);

	audio_tid = k_thread_create(&audio_thread_data, audio_stack_area, 
								K_THREAD_STACK_SIZEOF(audio_stack_area), 
								(k_thread_entry_t)ac_subsys, 
								dev, NULL, NULL, AUDIO_PRIORITY, 0, K_MSEC(0));
    k_thread_name_set(&audio_thread_data, "Audio RX");

	return ret;
}

void ac_deinit_thread(){
	k_thread_join(audio_tid, K_MSEC(200));
	k_thread_abort(audio_tid);
}

int main(void) {
	printk("[LOG] Starting application...\n");
	int ret;

	/* Set up PDM microphone */
	printk("[LOG] PDM setup\n");

	ret = device_is_ready(pdm_dev);
	if (ret < 0) {
		printk("Device Microphone is not ready\n");
		return ret;
	}

	ret = device_is_ready(pressure_dev);
	if (ret < 0) {
		printk("Device CF Sensor is not ready\n");
		return ret;
	}

	printk("[LOG] Starting Audio thread\n");
	ret = ac_init_thread(pdm_dev);

	/* Set up pressure */
	printk("[LOG] Starting Pressure thread\n");
	k_work_queue_init(&my_work_q);
	k_work_queue_start(&my_work_q, my_stack_area,
					K_THREAD_STACK_SIZEOF(my_stack_area), MY_PRIORITY,
					NULL);
	k_timer_start(&pressure_timer, K_MSEC(5), K_MSEC(5));

	// Main while loop
	while(1){
		k_msleep(1000);
	}

	// Shouldn't reach here
	ac_deinit_thread();
	printk("[LOG] Ending application...");

	return ret;
}