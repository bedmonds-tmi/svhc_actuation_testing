#include "pdm_microphone.h"

// -----------------------------------
// PDM microphone
// -----------------------------------
K_MEM_SLAB_DEFINE_STATIC(mem_slab, BLOCK_SIZE, BLOCK_COUNT, 4);
RING_BUF_ITEM_DECLARE_POW2(ring_buf, 15);

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
// Audio process management
// -----------------------------------
k_tid_t audio_tid = NULL;
struct k_thread audio_thread_data;
K_THREAD_STACK_DEFINE(audio_stack_area, AUDIO_STACK_SIZE);

/**
 * @brief Audio processing loop, could be used for all actuation detection processing
 */
static void audio_loop(const struct device* dev){
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
		ret = ring_buf_put(&ring_buf, buffer, size); // Returns number of bytes copied to ring buffer
		// if (ret < size){
		// 	printk("RING BUFFER FULL");
		// }
		k_mem_slab_free(&mem_slab, (void *)buffer);
	}

	ret = pdm_stop(dev);
	if(ret < 0) {
		return;
	}
}

int audio_init_record_thread(const struct device* dev){
	int ret;
	ret = pdm_init(dev);

	audio_tid = k_thread_create(&audio_thread_data, audio_stack_area, 
								K_THREAD_STACK_SIZEOF(audio_stack_area), 
								(k_thread_entry_t)audio_loop, 
								dev, NULL, NULL, AUDIO_PRIORITY, 0, K_MSEC(0));
    k_thread_name_set(&audio_thread_data, "Audio RX");

	return ret;
}

void audio_deinit_record_thread(){
	k_thread_join(audio_tid, K_MSEC(200));
	k_thread_abort(audio_tid);
}

// -----------------------------------
// USB Audio
// -----------------------------------
NET_BUF_POOL_FIXED_DEFINE(usb_audio_pool, CONFIG_USB_MAX_NUM_TRANSFERS, USB_BLOCK_SIZE, 0, NULL);
static void process_data(const struct device* dev) {
	int ret;

    if (ring_buf_size_get(&ring_buf) < USB_BLOCK_SIZE) {
        return;
    }

    // Allocate a new buffer inside the memory pool
    struct net_buf *buf = net_buf_alloc_fixed(&usb_audio_pool, K_NO_WAIT);

    if (buf != NULL) {
        // Copy data from ring buffer to temp buffer
        uint8_t tmp_buf[USB_BLOCK_SIZE];
        ring_buf_get(&ring_buf, tmp_buf, USB_BLOCK_SIZE);
		// Copy data from temp buffer to the buffer in the memory pool
        net_buf_add_mem(buf, tmp_buf, USB_BLOCK_SIZE);

        ret = usb_audio_send(dev, buf, USB_BLOCK_SIZE);
        if (ret < 0) {
            printk("[ERROR] USB Audio data_request_cb error");
            net_buf_unref(buf);
        }
    } else {
        printk("[ERROR] USB net_buf alloc error");
    }
}

static const struct usb_audio_ops mic_ops = {
	.data_request_cb = process_data
};

int usb_audio_init(const struct device* dev) {
	int ret;

	usb_audio_register(dev, &mic_ops);

	ret = usb_enable(NULL);
	if (ret != 0) {
        printk("[ERROR] USB Audio failed to enable");
		return -1;
	}

	return 0;
}
