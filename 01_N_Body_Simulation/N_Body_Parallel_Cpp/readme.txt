Compile:
g++ -I../../dependencies/ n_body.cpp -o n_body.exe

Run:
./n_body.exe --particles=1024 --mass=1.0e9 --dt=1.0 --iterations=20