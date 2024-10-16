## 此文件夹内final_main.cpp的目的是实现SDF场转stl，代替原来的matlab代码
C++17直接编译运行
argv[0]:cpp编译出来的可执行文件
argv[1]:csv转txt的点的数据，txt的路径
argv[2]:输出stl，stl的路径
e.g:
```
g++ -std=c++17 final_main.cpp -o final_main
./final_main /Users/merinomo/Documents/代码/SDFGen_project/Merino/dpp_txt.txt /Users/merinomo/Documents/代码/SDFGen_project/Merino/testdpp.stl
```