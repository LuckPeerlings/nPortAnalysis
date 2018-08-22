@echo off
curl.exe -s -k --user "api:key-1a232508ff0f409602a636381e6fcfd9" ^
 https://api.mailgun.net/v3/sandbox9abb08b487b2466ea1813f1da2005211.mailgun.org/messages ^
 -F from="The seven MWL dwarves <dwarves@sandbox9abb08b487b2466ea1813f1da2005211.mailgun.org>" ^
 -F to="jorgarti@gmail.com" ^
 -F to="luck@kth.se" ^
 -F subject="Measurement complete!" ^
 -F text="It is done, come down to the mine!"