@echo off
set /p scenario=Please enter the scenario number: 
echo ---------------------------------------------
set /p fiter=Please enter the first iteration number: 
echo ---------------------------------------------
set /p niter=Please enter the number of iterations: 
echo ---------------------------------------------
set /p mcmc=Please enter the MCMC duration (short, medium, or long):
echo ---------------------------------------------

echo Okay, I will now do your job for you.

"C:\Users\bstaton\Documents\R\R-4.0.2\bin\Rscript.exe" A-run-sims.R %scenario% %fiter% %niter% %mcmc%

pause
