#include <zephyr/kernel.h>
#include <zephyr/audio/dmic.h>
#include <zephyr/usb/usb_device.h>
#include <zephyr/usb/class/usb_audio.h>

#include <zephyr/sys/ring_buffer.h>

#define SAMPLE_FREQUENCY   16000
#define SAMPLE_BIT_WIDTH   16
#define BYTES_PER_SAMPLE   sizeof(int16_t)
#define NUMBER_OF_CHANNELS 2

#define SAMPLES_PER_BLOCK      ((SAMPLE_FREQUENCY / 50) * NUMBER_OF_CHANNELS)
#define SAMPLES_PER_BLOCK_OPUS ((SAMPLE_FREQUENCY / 50) * NUMBER_OF_CHANNELS)
#define INITIAL_BLOCKS         10
#define READ_TIMEOUT           1000

#define BLOCK_SIZE      (BYTES_PER_SAMPLE * SAMPLES_PER_BLOCK)
#define BLOCK_SIZE_OPUS (BYTES_PER_SAMPLE * SAMPLES_PER_BLOCK_OPUS)
#define BLOCK_COUNT     (INITIAL_BLOCKS)

// -----------------------------------
// Recording thread handling
// -----------------------------------
#define AUDIO_STACK_SIZE 5000
#define AUDIO_PRIORITY   -2
int audio_init_record_thread(const struct device* dev);
void audio_deinit_record_thread();

// -----------------------------------
// USB Audio handling, only needed for testing
// -----------------------------------
#define USB_BLOCK_SIZE 192 // Buffer size required for USB Audio
int usb_audio_init(const struct device *dev);