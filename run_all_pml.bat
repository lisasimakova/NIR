@echo off
setlocal
set PYMOL_DIR=C:\Users\Redmi\Desktop\NIR\mutation\file_PML
set LOGFILE=C:\Users\Redmi\Desktop\NIR\mutation\output_PML.txt

echo PyMOL Batch Run > "%LOGFILE%"
for %%F in ("%PYMOL_DIR%\*.pml") do (
    echo Running %%F >> "%LOGFILE%"
    pymol -cq "%%F" >> "%LOGFILE%" 2>&1
)
endlocal
pause
