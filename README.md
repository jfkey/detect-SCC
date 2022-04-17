# detect-SCC
Code for 

# Requirements

- Boost C++ Libraries (Version 1.71.0)
- Visual Studio 2017 (MSVC compiler)

# Usage


- Running our single incremental 
 

`-f  %dir%  -t "inc" -b 1 -a "tc_ahrsz" -n %nn%  -c 0  -p "comb"  -g "indpp_startC.txt"  -s true` 

`-f  %dir%  -t "batch" -b 1 -a "bat_cd4cg" -n %nn%  -c 0  -p "comb"  -g "indpp_startC.txt"  -s true` 

- Running batch single incremental 

`-f  %dir%  -t "batch" -b %batsize% -a "bat_cd4cg" -n %nn% -c  -1  -p "comb"  -g "indpp_startC.txt"  -s true`
