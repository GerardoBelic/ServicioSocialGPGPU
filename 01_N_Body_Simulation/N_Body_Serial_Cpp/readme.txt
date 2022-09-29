Compile:
g++ -I../../dependencies/include/ -std=c++20 -O2 -s n_body_serial.cpp -o n_body.exe

Run:
./n_body.exe --particles=1024 --mass=1.0e9 --dt=1.0 --iterations=200