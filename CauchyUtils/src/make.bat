R CMD SHLIB cauchyUtils.cpp -O2 cauchyUtils_main.cpp -o cauchyUtils.dll || exit
Rscript test.R 2>ff 1>out