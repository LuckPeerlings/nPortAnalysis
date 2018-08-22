% This function decodes error messages when the returned "status" 
% value from VT1432A module functions is not zero.
 
% Rev 1-8-98
 

function [session,status] = vt1432_check_status(session,status)
if(status)
   [status2,message]=vt1432('error_message',session,status);
   disp(['VISA Error: ' message]);
   [status2,details]=vt1432('errorDetails', session, 128);
   disp(['VISA Error: ' details]);
   disp(' ');
   disp('Paused... Hit CTRL-C to exit or Return to Continue');
   pause;
end










 