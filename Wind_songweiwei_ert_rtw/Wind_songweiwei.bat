set MATLAB=F:\MATLAB\R2014b

cd .

if "%1"=="" (F:\MATLAB\R2014b\bin\win64\gmake -f Wind_songweiwei.mk all) else (F:\MATLAB\R2014b\bin\win64\gmake -f Wind_songweiwei.mk %1)
@if errorlevel 1 goto error_exit

exit /B 0

:error_exit
echo The make command returned an error of %errorlevel%
An_error_occurred_during_the_call_to_make
