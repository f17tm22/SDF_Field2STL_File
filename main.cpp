/*
#include "ndgrid.cpp"
#ifndef BOWYERWATSON_CPP
#define BOWYERWATSON_CPP
#include "bowyerWatson.cpp"
#endif // BOWYERWATSON_CPP
*/
#include <cmath>
#include <algorithm>
#include <codecvt>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <time.h>
#include <unordered_set>
#include <vector>
#include <stack>
using std::pair;
using std::stack;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::istringstream;
using std::ofstream;
using std::string;
using std::stringstream;
using std::unordered_set;
using std::vector;
// according to the edge table, find the edges that intersect the
// isosurface
int edgeTable[256] = {
    0x0,   0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f,
    0xb06, 0xc0a, 0xd03, 0xe09, 0xf00, 0x190, 0x99,  0x393, 0x29a, 0x596, 0x49f,
    0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90, 0x230,
    0x339, 0x33,  0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936,
    0xe3a, 0xf33, 0xc39, 0xd30, 0x3a0, 0x2a9, 0x1a3, 0xaa,  0x7a6, 0x6af, 0x5a5,
    0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0, 0x460, 0x569,
    0x663, 0x76a, 0x66,  0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a,
    0x963, 0xa69, 0xb60, 0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff,  0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453,
    0x55a, 0x256, 0x35f, 0x55,  0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53,
    0x859, 0x950, 0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,  0xfcc,
    0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 0x8c0, 0x9c9, 0xac3, 0xbca,
    0xcc6, 0xdcf, 0xec5, 0xfcc, 0xcc,  0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9,
    0x7c0, 0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x55,
    0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650, 0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6,
    0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0xff,  0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f,
    0x66,  0x76a, 0x663, 0x569, 0x460, 0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af,
    0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa,  0x1a3, 0x2a9, 0x3a0, 0xd30,
    0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636,
    0x13a, 0x33,  0x339, 0x230, 0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895,
    0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99,  0x190, 0xf00, 0xe09,
    0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a,
    0x203, 0x109, 0x0};

struct sdfPoint {
  double x, y, z, sdfValue;
  // add default construct function
  sdfPoint(double x = 0.0f, double y = 0.0f, double z = 0.0f,
           double sdfValue = 0.0f)
      : x(x), y(y), z(z), sdfValue(sdfValue) {}
  bool operator==(const sdfPoint &other) const {
    return x == other.x && y == other.y && z == other.z &&
           sdfValue == other.sdfValue;
  }
};

struct Triangle {
  sdfPoint p1, p2, p3;
  Triangle(sdfPoint p1, sdfPoint p2, sdfPoint p3) : p1(p1), p2(p2), p3(p3) {}
  bool operator==(const Triangle &other) const {
    return p1 == other.p1 && p2 == other.p2 && p3 == other.p3;
  }
  bool sharesVertex(const Triangle &other) const {
    return p1 == other.p1 || p1 == other.p2 || p1 == other.p3 ||
           p2 == other.p1 || p2 == other.p2 || p2 == other.p3 ||
           p3 == other.p1 || p3 == other.p2 || p3 == other.p3;
  }
};

struct Edge {
  sdfPoint p1, p2;
  Edge(sdfPoint p1, sdfPoint p2) : p1(p1), p2(p2) {}
  bool operator==(const Edge &other) const {
    return p1 == other.p1 && p2 == other.p2;
  }
  bool sharesVertex(const Edge &other) const {
    return p1 == other.p1 || p1 == other.p2 || p2 == other.p1 || p2 == other.p2;
  }
};

struct Tetrahedron {
  sdfPoint p1, p2, p3, p4;
  Tetrahedron(sdfPoint p1, sdfPoint p2, sdfPoint p3, sdfPoint p4) : p1(p1), p2(p2), p3(p3), p4(p4) {}
  bool operator==(const Tetrahedron &other) const {
    return p1 == other.p1 && p2 == other.p2 && p3 == other.p3 && p4 == other.p4;
  }
};

vector<vector<double> > data;
vector<double> x;
vector<double> y;
vector<double> z;
vector<double> sdfValues;
vector<sdfPoint> points;
vector<double> uniqueX;
vector<double> uniqueY;
vector<double> uniqueZ;
vector<double> ux;
vector<double> uy;
vector<double> uz;

double vmin(const std::vector<double> &vec) {
  double vmin = 1e9;
  for (auto &value : vec) {
    if (vmin >= value)
      vmin = value;
  }
  return vmin;
}

double vmax(const std::vector<double> &vec) {
  double vmax = -1e9;
  for (auto &value : vec) {
    if (vmax <= value)
      vmax = value;
  }
  return vmax;
}

double fmin(double x, double y){ return x>y?y:x; }
double fmax(double x, double y){ return x<y?y:x; }

void readmatrix_txt(const string &txtFilePath, vector<vector<double> > &data) {
  // read data from csv file
  ifstream txt_data(txtFilePath, ios::in);

  if (!txt_data.is_open()) {
    cout << "Error: opening file fail" << endl;
    exit(1);
  } else {
    string line;
    vector<string> words;
    string word;
    vector<double> layer;
    txt_data.imbue(std::locale(txt_data.getloc()));
    // skip the first three lines at the beginning
    for (int i = 0; i < 3; ++i) {
      getline(txt_data, line);
    }

    // read data in rows
    while (getline(txt_data, line)) {
      istringstream sin(line);
      words.clear();
      // 将字符串流sin中的字符读到字符串数组words中，以逗号为分隔符
      vector<double> row; // 存储每一行的浮点数值
      stringstream ss(line);
      string value;

      while (getline(ss, value, ' ')) { // 使用逗号分隔每个字符串值
        try {
          if (value.find('e') != std::string::npos ||
              value.find('E') != std::string::npos) {
            double num =
                std::stod(value); // 使用std::stod转换科学计数法表示的值
            row.push_back(static_cast<double>(num));
          } else {
            double num =
                double(std::stof(value)); // 使用std::stof转换常规表示法的值
            row.push_back(num);
          }
        } catch (const std::invalid_argument &e) {
          cerr << "Error: Invalid argument encountered while converting "
                  "string to double. Invalid value: "
               << value << std::endl;
          continue;
        } catch (const std::out_of_range &e) {
          cerr << "Error: Out of range value encountered while "
                  "converting string to double."
               << std::endl;
          continue;
        }
      }

      data.push_back(row); // 将每一行的浮点数值存储到data中
      row.clear();         // 清空当前行的数据，准备读取下一行
    }
    // 打印转换后的浮点数值
    /*
    for (const auto &row : data) {
      for (const auto &num : row) {
        std::cout << num << " ";
      }
      std::cout << std::endl;
    }
    */
    txt_data.close();
  }
  cout << "...readmatrix_txt running complete." << endl;
}
void readmatrix_csv(const string &csvFilePath, vector<vector<double> > &data) {
  // read data from csv file
  ifstream csv_data(csvFilePath, ios::in);

  if (!csv_data.is_open()) {
    cout << "Error: opening file fail" << endl;
    exit(1);
  } else {
    string line;
    vector<string> words;
    string word;
    vector<double> layer;
    csv_data.imbue(std::locale(csv_data.getloc()));
    // skip the first three lines at the beginning
    for (int i = 0; i < 3; ++i) {
      getline(csv_data, line);
    }

    // read data in rows
    istringstream sin(line);
    words.clear();
    // 将字符串流sin中的字符读到字符串数组words中，以逗号为分隔符
    vector<double> row; // 存储每一行的浮点数值
    while (getline(sin, word, ',')) {
      stringstream ss(line);
      string value;

      while (std::getline(ss, value, ',')) { // 使用逗号分隔每个字符串值
        try {
          if (value.find('e') != std::string::npos ||
              value.find('E') != std::string::npos) {
            double num =
                std::stod(value); // 使用std::stod转换科学计数法表示的值
            row.push_back(static_cast<double>(num));
          } else {
            double num = std::stof(value); // 使用std::stof转换常规表示法的值
            row.push_back(num);
          }
        } catch (const std::invalid_argument &e) {
          cerr << "Error: Invalid argument encountered while converting "
                  "string to double. Invalid value: "
               << value << std::endl;
          continue;
        } catch (const std::out_of_range &e) {
          cerr << "Error: Out of range value encountered while "
                  "converting string to double."
               << std::endl;
          continue;
        }
      }

      data.push_back(row); // 将每一行的浮点数值存储到data中
      row.clear();         // 清空当前行的数据，准备读取下一行
    }
    // 打印转换后的浮点数值
    /*
    for (const auto &row : data) {
      for (const auto &num : row) {
        std::cout << num << " ";
      }
      std::cout << std::endl;
    }
    */
    csv_data.close();
  }
  cout << "...readmatrix_csv running complete." << endl;
}
// distinct whether the point is inside the circumcircle
bool isPointInsideCircumcircle(const sdfPoint& point, const Triangle& triangle) {
  // calculate the center of the circumcircle and radius
  float ax = triangle.p1.x;
  float ay = triangle.p1.y;
  float az = triangle.p1.z;
  float bx = triangle.p2.x;
  float by = triangle.p2.y;
  float bz = triangle.p2.z;
  float cx = triangle.p3.x;
  float cy = triangle.p3.y;
  float cz = triangle.p3.z;
  float d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
  float ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d;
  float uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d;
  float uz = ((ax * ax + ay * ay) * (bz - cz) + (bx * bx + by * by) * (cz - az) + (cx * cx + cy * cy) * (az - bz)) / d;
  float radius = sqrt((ax - ux) * (ax - ux) + (ay - uy) * (ay - uy) + (az - uz) * (az - uz));

  // distinct whether the point is inside the circumcircle
  float distance = sqrt((point.x - ux) * (point.x - ux) + (point.y - uy) * (point.y - uy) + (point.z - uz) * (point.z - uz));
  return distance <= radius;
}

bool isTriangleContainsPoint(const Triangle& triangle, const sdfPoint& point) {
  // use barycentric coordinates to distinct whether the point is inside the triangle
  float ax = triangle.p1.x;
  float ay = triangle.p1.y;
  float az = triangle.p1.z;
  float bx = triangle.p2.x;
  float by = triangle.p2.y;
  float bz = triangle.p2.z;
  float cx = triangle.p3.x;
  float cy = triangle.p3.y;
  float cz = triangle.p3.z;
  float dx = point.x;
  float dy = point.y;
  float dz = point.z;
  float d = (ay - cy) * (ax - cx) + (cx - bx) * (cy - by);
  float lambda1 = ((ay - cy) * (dx - cx) + (cx - bx) * (dy - cy)) / d;
  float lambda2 = ((cy - ay) * (dx - cx) + (ax - cx) * (dy - cy)) / d;
  float lambda3 = 1 - lambda1 - lambda2;
  return lambda1 >= 0 && lambda2 >= 0 && lambda3 >= 0;
}

// bowyerWatson algorithm O(nlogn)
vector<Tetrahedron> divideAndConquerTriangulation3D(const vector<sdfPoint>& points) {
  vector<Tetrahedron> result;
  stack<pair<vector<sdfPoint>, int> > stack;
  pair<vector<sdfPoint>, int> p(points, 0);
  stack.push(p);

  while (!stack.empty()) {
    auto [currentPoints, layer] = stack.top();
    stack.pop();

    //cout << "Layer: " << layer << std::endl;

    if (currentPoints.empty() || currentPoints.size() == 1) {
        // Base case: No triangulation needed
        continue;
    }

    if (currentPoints.size() == 2) {
        // Base case: Connect the two points to form an edge
        Tetrahedron tetrahedron(currentPoints[0], currentPoints[1], currentPoints[1], currentPoints[1]);
        // Process the tetrahedron...
        result.push_back(tetrahedron);
        continue;
    }

    // Divide the current points into two halves
    int mid = currentPoints.size() >> 1;
    vector<sdfPoint> leftPoints(currentPoints.begin(), currentPoints.begin() + mid);
    vector<sdfPoint> rightPoints(currentPoints.begin() + mid, currentPoints.end());

    // Push the left and right halves into the stack for further processing
    pair<vector<sdfPoint>, int> p1(leftPoints, layer+1);
    pair<vector<sdfPoint>, int> p2(rightPoints, layer+1);
    stack.push(p1);
    stack.push(p2);
    //cout << "Layer:" << layer << " finished" << endl;
  }
  return result;
}

void writeSTLASCII(const string& stlFilePath, const vector<Triangle>& triangles) {
    ofstream f(stlFilePath, ios::out | ios::trunc);
    if(!f.is_open()){
      cerr << "stlFile open failed: " << stlFilePath << endl;
      exit(1);
    } else {
      f.close();
      cout << "clear stlFile successfully: " << stlFilePath << endl;
    }
    ofstream file(stlFilePath, ios::out | ios::app);
    if (!file) {
        cerr << "Failed to open file for writing: " << stlFilePath << endl;
        exit(1);
    }
    
    // 写入头部信息
    file << "solid STLModel" << endl;
    
    // 逐个写入三角形面片的数据
    for (const auto& triangle : triangles) {
        // 写入法线向量
        file << "facet normal 0 0 0" << endl;
        file << "outer loop" << std::endl;
        
        // 写入三个顶点坐标
        file << "vertex " << triangle.p1.x << " " << triangle.p1.y << " " << triangle.p1.z << std::endl;
        file << "vertex " << triangle.p2.x << " " << triangle.p2.y << " " << triangle.p2.z << std::endl;
        file << "vertex " << triangle.p3.x << " " << triangle.p3.y << " " << triangle.p3.z << std::endl;
        
        file << "endloop" << endl;
        file << "endfacet" << endl;
    }
    
    file << "endsolid STLModel" << endl;
    
    file.close();
}

vector<sdfPoint> marchingCubes() {
  // create 3D grid
  double enlarge_rate = 1.000000;
  int gridSizeX = ux.size() * enlarge_rate;
  int gridSizeY = uy.size() * enlarge_rate;
  int gridSizeZ = uz.size() * enlarge_rate;
  cout << "grid size: " << gridSizeX << ' ' << gridSizeY << ' ' << gridSizeZ << endl;
  std::vector<std::vector<std::vector<double> > > sdfGrid(
      gridSizeX, std::vector<std::vector<double> >(
                     gridSizeY, std::vector<double>(gridSizeZ)));

  // fill the SDF value into the grid

  int xIndex, yIndex, zIndex;
  //int xidx, yidx, zidx;
  // return;
  //cout << "breakpoint1" << endl;
  //int cnt = 0;
  /*
  ofstream f("/Users/merinomo/Documents/代码/SDFGen_project/Merino/test.txt",
  ofstream::out | ofstream::trunc); if (f.is_open()) {
    // clear file
    f.close();
    cout << "...clear f successfully"<< endl;
  } else {
    cout << "...Fail to open f"<< endl;
  }
  ofstream ff("/Users/merinomo/Documents/代码/SDFGen_project/Merino/test.txt",
  ofstream::out); for (const auto &point : points) {
    ++cnt;
    ff << cnt << ':' << point.x << ' ' << point.y << ' ' << point.z << ' ' <<
  point.sdfValue << endl;
  }
  */
  //cnt = 0;
  for (const auto &point : points) {
    //++cnt;
    // cout << point.x << ' ' << point.y << ' ' << point.z << ' ' <<
    // point.sdfValue << endl; cout << *ux.begin() << ' ' << *ux.end() << ' ' <<
    // *lower_bound(ux.begin(), ux.end(), point.x) << endl; cout << *uy.begin()
    // << ' ' << *uy.end() << ' ' << *lower_bound(uy.begin(), uy.end(), point.y)
    // << endl; cout << *uz.begin() << ' ' << *uz.end() << ' ' <<
    // *lower_bound(uz.begin(), uz.end(), point.z) << endl;
    xIndex = std::distance(ux.begin(),
                           std::lower_bound(ux.begin(), ux.end(), point.x));
    yIndex = std::distance(uy.begin(),
                           std::lower_bound(uy.begin(), uy.end(), point.y));
    zIndex = std::distance(uz.begin(),
                           std::lower_bound(uz.begin(), uz.end(), point.z));
    sdfGrid[xIndex][yIndex][zIndex] = point.sdfValue;
    // cout << cnt << ':' << xIndex << ' ' << yIndex << ' ' << zIndex << ' ' <<
    // sdfGrid[xIndex][yIndex][zIndex] << endl;
  }
  //cout << "breakpoint2" << endl;
  /*
  for (int i = 0; i < xIndex; ++i) {
    for(int j = 0; j < yIndex; ++j) {
      for(int k = 0; k < zIndex; ++k) {
        cout << i << ',' << j << ',' << k << ':' << sdfGrid[i][j][k] << endl;
      }
      cout << endl;
    }
    cout << endl;
  }
  */
  cout << "...create 3D grid complete." << endl;
  // marchingCubes algorithm
  // search each cell in the grid
  vector<sdfPoint> vertices;
  for (int x = 0; x < gridSizeX - 1; ++x) {
    for (int y = 0; y < gridSizeY - 1; ++y) {
      for (int z = 0; z < gridSizeZ - 1; ++z) {
        // get the coordinates and SDF values of the 8 vertices of the current
        // cell

        double sdfValues[8] = {sdfGrid[x][y][z],
                               sdfGrid[x + 1][y][z],
                               sdfGrid[x + 1][y + 1][z],
                               sdfGrid[x][y + 1][z],
                               sdfGrid[x][y][z + 1],
                               sdfGrid[x + 1][y][z + 1],
                               sdfGrid[x + 1][y + 1][z + 1],
                               sdfGrid[x][y + 1][z + 1]};
        // distinct whether the cell is intersecting the surface
        int cubeIndex = 0;
        if (sdfValues[0] < 0)
          cubeIndex |= 1;
        if (sdfValues[1] < 0)
          cubeIndex |= 2;
        if (sdfValues[2] < 0)
          cubeIndex |= 4;
        if (sdfValues[3] < 0)
          cubeIndex |= 8;
        if (sdfValues[4] < 0)
          cubeIndex |= 16;
        if (sdfValues[5] < 0)
          cubeIndex |= 32;
        if (sdfValues[6] < 0)
          cubeIndex |= 64;
        if (sdfValues[7] < 0)
          cubeIndex |= 128;
        // according to the edge table, get the edges that intersect the
        // isosurface in the current cell
        int edgeIndices[12] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
        if (edgeTable[cubeIndex] & 1)
          edgeIndices[0] = 0;
        if (edgeTable[cubeIndex] & 2)
          edgeIndices[1] = 1;
        if (edgeTable[cubeIndex] & 4)
          edgeIndices[2] = 2;
        if (edgeTable[cubeIndex] & 8)
          edgeIndices[3] = 3;
        if (edgeTable[cubeIndex] & 16)
          edgeIndices[4] = 4;
        if (edgeTable[cubeIndex] & 32)
          edgeIndices[5] = 5;
        if (edgeTable[cubeIndex] & 64)
          edgeIndices[6] = 6;
        if (edgeTable[cubeIndex] & 128)
          edgeIndices[7] = 7;
        if (edgeTable[cubeIndex] & 256)
          edgeIndices[8] = 8;
        if (edgeTable[cubeIndex] & 512)
          edgeIndices[9] = 9;
        if (edgeTable[cubeIndex] & 1024)
          edgeIndices[10] = 10;
        if (edgeTable[cubeIndex] & 2048)
          edgeIndices[11] = 11;

        // according to the index of edges and the interpolation, calculate the
        // coordinates of the vertice
        // There exists problem...
        for (int i = 0; i < 12; ++i) {
          if (edgeIndices[i] != -1 && edgeIndices[(i + 1) % 12] != -1) {
            if (sdfValues[edgeIndices[i]] != 0 &&
                sdfValues[edgeIndices[(i + 1) % 12]] != 0) {
              double t = 1.000000*(sdfValues[edgeIndices[i]] /
                         (sdfValues[edgeIndices[i]] -
                          sdfValues[edgeIndices[(i + 1) % 12]]));
              //cout << "t:" << t << endl;
              if(t > 1 || t < -1 || t == 0){
                continue;
              }
              sdfPoint newVertex(0, 0, 0, 0);
              newVertex.x = 1.000000*x + t;
              newVertex.y = 1.000000*y + t;
              newVertex.z = 1.000000*z + t;
              newVertex.sdfValue = 1.000000*sdfValues[0] + t;
              //cout << newVertex.x << ' ' << newVertex.y << ' ' << newVertex.z << ' ' << newVertex.sdfValue << endl;
              if (std::isinf(newVertex.x) || std::isinf(newVertex.y) || 
                  std::isinf(newVertex.z) || !newVertex.x || 
                  !newVertex.y || !newVertex.z || std::isnan(newVertex.x) || 
                  std::isnan(newVertex.y) || std::isnan(newVertex.z)) {
                // cout << edgeIndices[i] << ' ' << cubeIndex << ' ' <<
                // edgeTable[cubeIndex] << ' ' << endl;
                continue;
              }
              vertices.push_back(newVertex);
              //cout << newVertex.x << ' ' << newVertex.y << ' ' << newVertex.z << ' ' << newVertex.sdfValue << endl;
            }
          }
        }
        //cout << sdfValues[0] << ' ';
      }
    }
    // loading...
    /*
    int process = ceil(50.00 * (x + 1) / (gridSizeX - 1));
    for (int i = 0; i < process; ++i) {
      cout << '*';
    }
    for (int i = 0; i < 50 - process; ++i) {
      cout << ' ';
    }
    cout << (100.00 * (x + 1) / (gridSizeX - 1)) << '%' << endl;
    */
  }
  // print the unique coordinate
  /*
  for(const Point& p : points){
    cout << p.x << " " << p.y << " " << p.z << endl;
  }
  for (const auto& num : uniqueX) {
    cout << num << " ";
  }
  cout << endl;
  for (const auto& num : uniqueY) {
    cout << num << " ";
  }
  cout << endl;
  for (const auto& num : uniqueZ) {
    cout << num << " ";
  }
  cout << endl;
  */
  cout << "...marchingCubes running complete." << endl;
  return vertices;
}

// extract the unique value from vectors x, y, z
vector<double> unique(const vector<double> &vec) {
  std::unordered_set<double> seen;
  vector<double> uniqueVec;
  for (const auto &value : vec) {
    if (seen.find(value) == seen.end()) {
      seen.insert(value);
      uniqueVec.push_back(value);
    }
  }
  return uniqueVec;
}

std::string getFileExtension(const std::string &filePath) {
  size_t dotIndex = filePath.find_last_of(".");
  if (dotIndex != std::string::npos && dotIndex < filePath.length() - 1) {
    return filePath.substr(dotIndex + 1);
  }
  return "";
}

int main(int argc, char *argv[]) {
  clock_t start, end;
  start = clock();
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " <input file>" << endl;
  }
  ofstream file(argv[2], ofstream::out | ofstream::trunc);
  if (file.is_open()) {
    // clear file
    file.close();
    cout << "...clear file successfully" << endl;
  } else {
    cout << "...Fail to open file" << endl;
  }
  if (getFileExtension(argv[1]) == "txt") {
    readmatrix_txt(argv[1], data);
  } else if (getFileExtension(argv[1]) == "csv") {
    readmatrix_csv(argv[1], data);
  }
  // cout data in order to test whether data is correctly read in
  /*
  for (const auto& layer : data) {
      for (const auto& value : layer) {
          cout << value << " ";
      }
      cout << endl;
  }
  */
  for (const auto &row : data) {
    x.push_back(row[0]);
    y.push_back(row[1]);
    z.push_back(row[2]);
    sdfValues.push_back(row[3]);
    // cout << row[0] << ' ' << row[1] << ' ' << row[2] << ' ' << row[3] <<
    // endl;
    sdfPoint p(row[0], row[1], row[2], row[3]);
    // cout << p.x << " " << p.y << " " << p.z << " " << p.sdfValue << endl;
    points.push_back(p);
  }
  /* points correct.
  for(auto& point : points){
    cout << point.x << " " << point.y << " " << point.z << " " << point.sdfValue
  << endl;
  }
  */
  cout << "...extract data complete." << endl;
  // find unique coordinate from x, y, z
  uniqueX = unique(x);
  uniqueY = unique(y);
  uniqueZ = unique(z);
  // remove the coordinate 5e8
  for (const auto &x : uniqueX) {
    if (x >= 4e8)
      continue;
    ux.push_back(x);
    // cout << x << endl;
  }
  for (const auto &y : uniqueY) {
    if (y >= 4e8)
      continue;
    uy.push_back(y);
    // cout << y << endl;
  }
  for (const auto &z : uniqueZ) {
    if (z >= 4e8)
      continue;
    uz.push_back(z);
    // cout << z << endl;
  }

  // cout the unique coordinate
  /*
  cout << "Unique X values: ";
  for (double val : ux) std::cout << val << " ";
  cout << std::endl;

  cout << "Unique Y values: ";
  for (double val : uy) std::cout << val << " ";
  cout << std::endl;

  cout << "Unique Z values: ";
  for (double val : uz) std::cout << val << " ";
  cout << std::endl;
  */

  // get the range of coordinate of x, y, z
  double xMin = vmin(ux);
  double xMax = vmax(ux);
  double yMin = vmin(uy);
  double yMax = vmax(uy);
  double zMin = vmin(uz);
  double zMax = vmax(uz);
  cout << "min: " << xMin << " " << yMin << " " << zMin << endl;
  cout << "max: " << xMax << " " << yMax << " " << zMax << endl;
  // calculate the center points of x, y, z axies
  double centerX = (xMax + xMin) / 2;
  double centerY = (yMax + yMin) / 2;
  double centerZ = (zMax + zMin) / 2;
  cout << "center: " << centerX << " " << centerY << " " << centerZ << endl;
  // translate the coordinates and ensure that model is placed at the center of
  // grid
  for (size_t i = 0; i < ux.size(); ++i) {
    ux[i] -= centerX;
  }
  for (size_t i = 0; i < uy.size(); ++i) {
    uy[i] -= centerY;
  }
  for (size_t i = 0; i < uz.size(); ++i) {
    uz[i] -= centerZ;
  }
  // create 3D grid
  // array<vector<double>, 3> coordniate = ndgrid(uniqueX, uniqueY, uniqueZ);
  vector<sdfPoint> vertices = marchingCubes();
  /*
  for(auto& v: vertices){
    cout << v.x << ' ' << v.y << ' ' << v.z << ' ' << v.sdfValue << endl;
  }
  */
  vector<Tetrahedron> tetrahedron = divideAndConquerTriangulation3D(vertices);
  int cnt = 0;
  vector<Triangle> triangles;
  //transform tetrahedron into 
  for(auto& t : tetrahedron){
    Triangle t1(t.p1, t.p2, t.p3);
    Triangle t2(t.p1, t.p2, t.p4);
    Triangle t3(t.p1, t.p3, t.p4);
    Triangle t4(t.p2, t.p3, t.p4);
    triangles.push_back(t1);
    triangles.push_back(t2);
    triangles.push_back(t3);
    triangles.push_back(t4);
    /*
    cout << "tetrahedron" << cnt << ':' << endl;
    ++cnt;
    cout << "t1:" <<t.p1.x << ' ' << t.p1.y << ' ' << t.p1.z << ' ' << t.p1.sdfValue << endl;
    cout << "t2:" <<t.p2.x << ' ' << t.p2.y << ' ' << t.p2.z << ' ' << t.p2.sdfValue << endl;
    cout << "t3:" <<t.p3.x << ' ' << t.p3.y << ' ' << t.p3.z << ' ' << t.p3.sdfValue << endl;
    cout << "t4:" <<t.p4.x << ' ' << t.p4.y << ' ' << t.p4.z << ' ' << t.p4.sdfValue << endl;
    */
  }
  writeSTLASCII(argv[2], triangles);
  cout << "...final_main running complete." << endl;
  end = clock(); // 结束时间
  cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s"
       << endl; // 输出时间（单位：ｓ）
  return 0;
}
