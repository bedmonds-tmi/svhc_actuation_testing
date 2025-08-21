#pragma once
#include <zephyr/audio/dmic.h>

#define SAMPLE_FREQUENCY   16000
#define SAMPLE_BIT_WIDTH   16
#define BYTES_PER_SAMPLE   sizeof(int16_t)

#define SAMPLES_PER_BLOCK      2048
#define SAMPLES_PER_BLOCK_OPUS 2048
#define INITIAL_BLOCKS         10
#define READ_TIMEOUT           1000

#define BLOCK_SIZE      (BYTES_PER_SAMPLE * SAMPLES_PER_BLOCK)
#define BLOCK_SIZE_OPUS (BYTES_PER_SAMPLE * SAMPLES_PER_BLOCK_OPUS)
#define BLOCK_COUNT     (INITIAL_BLOCKS)