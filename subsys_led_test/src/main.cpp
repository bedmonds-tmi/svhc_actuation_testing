#include <zephyr/kernel.h>
#include <zephyr/logging/log.h>
#include <tmi/sub/led/led.hpp>

LOG_MODULE_REGISTER(app);

int main(void)
{
	using ledCmd = ssLed::Command::Type;
	subLed.init();
  subLed.start();
  subLed.cmd({ledCmd::SET_STATE, ssLed::State::BLINK});
  while(1) {
    k_msleep(1000);
  }
  return 0;
}