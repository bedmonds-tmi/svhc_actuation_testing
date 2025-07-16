#include <zephyr/kernel.h>
#include <zephyr/audio/dmic.h>
#include <zephyr/sys/ring_buffer.h>
#include <nrfx_clock.h>

#include "audio_settings.h"

/* ---------- PDM Microphone ---------- 
- mem_slab: a memory slab of BLOCK_COUNT blocks of size BLOCK_SIZE,
-- passed to the PDM module (via the dmic zephyr API) and filled via DMA
-- read one block at a time
*/
K_MEM_SLAB_DEFINE_STATIC(mem_slab, BLOCK_SIZE, BLOCK_COUNT, 4);
static const struct device *const dmic_dev = DEVICE_DT_GET(DT_NODELABEL(dmic_dev));

static inline int pdm_init()
{
  const struct device *dev
	int ret;

	if (!device_is_ready(dmic_dev)) {
		printk("%s is not ready", dmic_dev->name);
		return -1;
	}

	struct pcm_stream_cfg stream = {
		.pcm_rate = SAMPLE_FREQUENCY,
		.pcm_width = SAMPLE_BIT_WIDTH,
		.block_size = BLOCK_SIZE,
		.mem_slab  = &mem_slab,
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
			.max_pdm_clk_dc   = 60,
		},
		.streams = &stream,
		.channel = {
			.req_num_chan = 2,
			.req_num_streams = 1,
		},
	};

	cfg.channel.req_chan_map_lo = dmic_build_channel_map(1, 0, PDM_CHAN_LEFT) | dmic_build_channel_map(0, 0, PDM_CHAN_RIGHT);

	printk("PCM output rate: %u, channels: %u", cfg.streams[0].pcm_rate, cfg.channel.req_num_chan);

	ret = dmic_configure(dmic_dev, &cfg);
	if (ret < 0) {
		printk("Failed to configure the driver: %d", ret);
		return ret;
	}

	ret = dmic_trigger(dmic_dev, DMIC_TRIGGER_START);
	if (ret < 0) {
		printk("START trigger failed: %d", ret);
		return ret;
	}

	return 0;
}

/* ---------- Processing ---------- 
- ring_buf: a ring buffer that PDM blocks are loaded into
*/
RING_BUF_ITEM_DECLARE_POW2(ring_buf, 15);

int main(void)
{
	printk("PDM to USB Audio\n");
	int ret;

	/* Set up PDM microphone */
	printk("PDM setup\n");
	ret = pdm_init();
	if(ret < 0) {
		return ret;
	}
	printk("PDM setup complete\n");

	while(1){
		/* Read data from PDM microphone */
		uint8_t *buffer = NULL;
		uint32_t size = 0;

		ret = dmic_read(dmic_dev, 0, (void *)&buffer, &size, READ_TIMEOUT);
		if (ret < 0) {
			printk("read failed: %d\n", ret);
			return ret;
		}

		/* Copy data to ring buffer */
		ret = ring_buf_put(&usb_ring_buf, buffer, size); // Returns number of bytes copied to ring buffer
		k_mem_slab_free(&mem_slab, (void *)buffer);
	}

	ret = dmic_trigger(dmic_dev, DMIC_TRIGGER_STOP);
	if (ret < 0) {
		printk("STOP trigger failed: %d", ret);
		return ret;
	}

	return ret;
}