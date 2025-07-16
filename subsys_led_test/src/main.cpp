#include <zephyr/kernel.h>
#include <zephyr/logging/log.h>
#include <tmi/sub/led/led.hpp>

LOG_MODULE_REGISTER(app)

int main(void)
{
	subLed.init();
  subLed.start();
  subLed.blink();
  while(1) {
    k_msleep(1000);
  }
}