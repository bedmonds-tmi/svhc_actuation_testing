#include <zephyr/kernel.h>
#include <zephyr/timing/timing.h>

#include "pdm_microphone.h"
#include <tmi/drv/xgzp68/xgzp68.h>

static const struct device *const pdm_dev = DEVICE_DT_GET(DT_NODELABEL(pdm_dev));
static const struct device *const pressure_dev = DEVICE_DT_GET_ONE(vnd_xgzp68);
static const struct device *const usb_dev = DEVICE_DT_GET_ONE(usb_audio_mic);

void print_pressure() {
	float pressure;
	// int64_t us = k_uptime_ticks() * 1000000LL / CONFIG_SYS_CLOCK_TICKS_PER_SEC;
	pressure_get(pressure_dev, &pressure);
	// printk("%lld us, %f\n", us, (double)pressure);
	printk("prs%fprs\n", (double)pressure);
}

K_WORK_DEFINE(pressure_work, print_pressure);

void pressure_timer_handler(struct k_timer *dummy) {
    k_work_submit(&pressure_work);
}

K_TIMER_DEFINE(pressure_timer, pressure_timer_handler, NULL);

int main(void) {
	printk("[LOG] Starting application...\n");
	int ret;

	/* Set up PDM microphone */
	printk("[LOG] PDM setup\n");

	ret = device_is_ready(pdm_dev);
	if (ret < 0) {
		printk("Device USB Microphone is not ready\n");
		return ret;
	}

	printk("[LOG] Starting Audio thread\n");
	ret = audio_init_record_thread(pdm_dev);

	/* Set up USB Audio */
	printk("[LOG] USB Audio setup\n");

	ret = device_is_ready(usb_dev);
	if (ret < 0) {
		printk("Device USB Microphone is not ready\n");
		return ret;
	}
	ret = usb_audio_init(usb_dev);
	if(ret < 0) {
		return ret;
	}
	printk("[LOG] USB Audio setup complete\n");

	/* Set up pressure */
	pressure_start(pressure_dev);
	k_timer_start(&pressure_timer, K_MSEC(5), K_MSEC(5));

	// Main while loop
	while(1){
		k_msleep(1000);
	}

	// Shouldn't reach here
	audio_deinit_record_thread();
	printk("[LOG] Ending application...");

	return ret;
}