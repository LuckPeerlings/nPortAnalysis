function [status,answer] = mailgun(email)
%% Sends an email through Mailgun API and returns status and server answer
% Example of the expected email structure:
% email.APIkey  = 'key-1a232508ff0f409602a636381e6fcfd9';
% email.domain  = 'sandbox9abb08b487b2466ea1813f1da2005211.mailgun.org';
% email.from    = 'The seven MWL dwarves <dwarves@sandbox9abb08b487b2466ea1813f1da2005211.mailgun.org>';
% email.to      = 'jorgarti@gmail.com';
% email.subject = 'Measurements finished';
% email.text    = 'Come down to the mine, we are done!';

%% HTTPS request
old = cd([fileparts(mfilename('fullpath')) '/curl']);
[status,answer] = system(['curl -k -s --user "api:' email.APIkey...
                          '" https://api.mailgun.net/v3/' email.domain...
                          '/messages -F from="' email.from...
                          '" -F to="' email.to...
                          '" -F subject="' email.subject...
                          '" -F text="' email.text '"']);
cd(old)