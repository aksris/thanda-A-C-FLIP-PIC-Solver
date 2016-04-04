#include "macgriddata.h"
#include <iostream>

float CellSize = 1.f;
int gDimension[3] = {5, 5, 5};
float LERP(float val1, float val2, float t){
    return (1 - t) * val1 + t * val2;
}

MACGridData::MACGridData():
containerBounds(5.0f, 5.f, 5.f)
{
}

MACGridData::MACGridData(const MACGridData &val){

    containerBounds = val.containerBounds;
}

MACGridData& MACGridData::operator =(const MACGridData& val){
    if (this == &val){
        return *this;
    }

    data = val.data;
    containerBounds = val.containerBounds;
    return *this;
}

int MACGridData::getCellIndex(int i, int j, int k){

    int x = i;
    int y = j * gDimension[0];
    int z = k * gDimension[0] * gDimension[1];

    return x+y+z;
}

void MACGridData::MACGridDataInitialize(){
    containerBounds[0] = CellSize * gDimension[0];
    containerBounds[1] = CellSize*gDimension[1];
    containerBounds[2] = CellSize*gDimension[2];
    data.resize(gDimension[0]*gDimension[1]*gDimension[2]);
    mData.resize(gDimension[0]*gDimension[1]*gDimension[2]);
    std::fill(data.begin(), data.end(), 0.f);
    std::fill(mData.begin(), mData.end(), 0.f);
}

float& MACGridData::operator ()(int i, int j, int k){

    float ret = 0.f;
    if (i < 0 || j < 0 || k < 0 ||
            i > gDimension[0]-1 ||
            j > gDimension[1]-1 ||
            k > gDimension[2]-1) return ret;

    int x = i;
    int y = j * gDimension[0];
    int z = k * gDimension[0] * gDimension[1];

    return data[x+y+z];
}

void MACGridData::getCell(const vec3 &pt, int &i, int &j, int &k){
    vec3 posL = worldToLocal(pt);
//    vec3 minBounds = -containerBounds.x / 2.f,
    i = (int) (posL[0]/CellSize);
    j = (int) (posL[1]/CellSize);
    k = (int) (posL[2]/CellSize);
}

int MACGridData::getCellMark(int i, int j, int k){
    float ret = 0.f;
    if (i < 0 || j < 0 || k < 0 ||
            i > gDimension[0]-1 ||
            j > gDimension[1]-1 ||
            k > gDimension[2]-1) return ret;

    int x = i;
    int y = j * gDimension[0];
    int z = k * gDimension[0] * gDimension[1];

    return mData[x+y+z];
}

void MACGridData::setCell(int &i, int &j, int &k, const float val){
    int x = i;
    int y = j * gDimension[0];
    int z = k * gDimension[0] * gDimension[1];
    data[x+y+z] = val;
}

void MACGridData::setCellMark(int &i, int &j, int &k, const int val, bool mark){
    int x = i;
    int y = j * gDimension[0];
    int z = k * gDimension[0] * gDimension[1];
    mData[x+y+z] = val;
}

vec3 MACGridData::worldToLocal(const vec3& pt) const
{
   vec3 ret;
   ret[0] = min(max(0.0f, pt[0] - CellSize*0.5f), containerBounds[0]); //staggering
   ret[1] = min(max(0.0f, pt[1] - CellSize*0.5f), containerBounds[1]);
   ret[2] = min(max(0.0f, pt[2] - CellSize*0.5f), containerBounds[2]);
   return ret;
}

float MACGridData::interpolate(const vec3& pt)
{
    
//#define DEBUG

    vec3 pos = worldToLocal(pt);

    int i = (int) (pos[0]/CellSize);
    int j = (int) (pos[1]/CellSize);
    int k = (int) (pos[2]/CellSize);

    //   float scale = 1.0 / CellSize;
    float fract_partx = (pos[0] - i*CellSize);
    float fract_party = (pos[1] - j*CellSize);
    float fract_partz = (pos[2] - k*CellSize);
#ifdef DEBUG
    std::cout << "[thanda] InterpolateVelocity"<< std::endl;
#endif
    // questions -
    // 1: Shouldn't all the neighbors have extrapolated values for correct interpolation - Yes!
    //    This function is trusting ExtrapolateVelocity to make sure there are no non-existant terms.
    //    Or are we supposed to check here regardless?
    
    // 2: Are boundary cells supposed to be type SOLID? If so why, is cell 0 2 1 not being extrapolated to correctly?
    
    float v000 = (*this)(i >= gDimension[0]? gDimension[0] - 1 : i, j >= gDimension[1]? gDimension[1]    -1: j, k >= gDimension[2]? gDimension[2] -1 : k);
    float v010 = (*this)(i >= gDimension[0]? gDimension[0] - 1 : i, j+1 >= gDimension[1]? gDimension[1]  -1: j + 1,k >= gDimension[2]? gDimension[2]-1 : k);
    float lerp1 = LERP(v000, v010, fract_party);
#ifdef DEBUG
    std::cout << "[thanda] at idx " << i << " "<< j << " "<< k << " : "<< (*this)(i,j,k) << std::endl;
    std::cout << "[thanda] at idx " << i << " "<< j+1 << " "<< k << " : "<< (*this)(i,j+1,k) << std::endl;
    std::cout << "[thanda] lerp1 " << lerp1 << std::endl;
#endif

    float v100 = (*this)(i+1 >= gDimension[0]? gDimension[0] -1 : i+1, j >= gDimension[1]? gDimension[1] -1: j, k >= gDimension[2]? gDimension[2] -1: k);
    float v110 = (*this)(i+1 >= gDimension[0]? gDimension[0] -1 : i+1,j+1 >= gDimension[1]? gDimension[1]-1: j + 1,k >= gDimension[2]? gDimension[2] -1: k);
    float lerp2 = LERP(v100, v110, fract_party);
#ifdef DEBUG
    std::cout << "[thanda] at idx " << i+1 << " "<< j << " "<< k << " : "<< (*this)(i+1,j,k) << std::endl;
    std::cout << "[thanda] at idx " << i+1 << " "<< j+1 << " "<< k << " : "<< (*this)(i+1,j+1,k) << std::endl;
    std::cout << "[thanda] lerp2 " << lerp2 << std::endl;
#endif
    
    float v001 = (*this)(i >= gDimension[0]? gDimension[0] -1 : i, j >= gDimension[1]? gDimension[1]     -1: j, k+1 >= gDimension[2] ? gDimension[2] -1: k + 1);
    float v011 = (*this)(i >= gDimension[0]? gDimension[0] -1 : i,j+1 >= gDimension[1]? gDimension[1]-1: j + 1,k+1 >= gDimension[2] ? gDimension[2] -1: k+1);
    float lerp3 = LERP(v001, v011, fract_party);
#ifdef DEBUG
    std::cout << "[thanda] at idx " << i << " "<< j << " "<< k+1<< " : "<< (*this)(i,j,k+1) << std::endl;
    std::cout << "[thanda] at idx " << i << " "<< j+1 << " "<< k+1 << " : "<< (*this)(i,j+1,k+1) << std::endl;
    std::cout << "[thanda] lerp3 " << lerp3 << std::endl;
#endif
    
    float v101 = (*this)(i+1 >= gDimension[0]? gDimension[0] -1 : i+1,j >= gDimension[1]? gDimension[1]    -1: j,k+1 >= gDimension[2] ? gDimension[2] -1: k + 1);
    float v111 = (*this)(i+1 >= gDimension[0]? gDimension[0] -1 : i+1,j+1 >= gDimension[1]? gDimension[1]  -1: j + 1,k+1 >= gDimension[2] ? gDimension[2]-1 : k + 1);
    float lerp4 = LERP(v101, v111, fract_party);
#ifdef DEBUG
    std::cout << "[thanda] at idx " << i+1 << " "<< j << " "<< k+1 << " : "<< (*this)(i+1,j,k+1) << std::endl;
    std::cout << "[thanda] at idx " << i+1 << " "<< j+1 << " "<< k+1 << " : "<< (*this)(i+1,j+1,k+1) << std::endl;
    std::cout << "[thanda] lerp4 " << lerp4 << std::endl;
#endif

    float lerp5 = LERP (lerp1, lerp2, fract_partx);
    float lerp6 = LERP (lerp3, lerp4, fract_partx);
    float ret = LERP(lerp5, lerp6, fract_partz);
#ifdef DEBUG
    std::cout << "[thanda] lerp5 " << lerp5 << std::endl;
    std::cout << "[thanda] lerp6 " << lerp6 << std::endl;
    std::cout << "[thanda] interpolated value " << ret << std::endl;
#endif


    return ret;
}

MACGridDataX::MACGridDataX()
{
}

void MACGridDataX::MACGridDataInitialize(){
//    MACGridData::MACGridDataInitialize();
    containerBounds[0] = CellSize * (gDimension[0] + 1);
    containerBounds[1] = CellSize * gDimension[1];
    containerBounds[2] = CellSize * gDimension[2];
    data.resize((gDimension[0] + 1) * gDimension[1] * gDimension[2]);
    std::fill(data.begin(), data.end(), 0.f);
}

float& MACGridDataX::operator ()(int i, int j, int k){

    float ret = 0.f;
    if (i < 0 || i > gDimension[0]) return ret;
    if (j < 0) {
        j = 0;
    }
    if (k < 0){
        k = 0;
    }
    if (j > gDimension[1]){
        j = gDimension[1] - 1;
    }
    if (k > gDimension[2]){
        k = gDimension[2] - 1;
    }

    int x = i;
    int y = j * (gDimension[0] + 1);
    int z = k * (gDimension[0] + 1) * gDimension[1];

    return data[x+y+z];
}

vec3 MACGridDataX::worldToLocal(const vec3 &pt) const {
    vec3 ret;
    ret[0] = min(max(0.0f, pt[0]), containerBounds[0]); //staggering // its already at 0 face of x axis
    ret[1] = min(max(0.0f, pt[1] - CellSize*0.5f), containerBounds[1]);
    ret[2] = min(max(0.0f, pt[2] - CellSize*0.5f), containerBounds[2]);
    return ret;
}

void MACGridDataX::setCell(int &i, int &j, int &k, const float val){
    int x = i;
    int y = j * (gDimension[0] + 1);
    int z = k * (gDimension[0] + 1) * gDimension[1];
    data[x+y+z] = val;
}

void MACGridDataX::setCellAdd(int &i, int &j, int &k, const float val){
    int x = i;
    int y = j * (gDimension[0] + 1);
    int z = k * (gDimension[0] + 1) * gDimension[1];
    data[x+y+z] += val;
}

MACGridDataY::MACGridDataY()
{
}

void MACGridDataY::MACGridDataInitialize(){
    containerBounds[0] = CellSize * (gDimension[0]);
    containerBounds[1] = CellSize * (gDimension[1] + 1);
    containerBounds[2] = CellSize * gDimension[2];
    data.resize((gDimension[0]) * (gDimension[1] + 1) * gDimension[2]);
    std::fill(data.begin(), data.end(), 0.f);
}

float& MACGridDataY::operator ()(int i, int j, int k){

    float ret = 0.f;
    if (j < 0 || j > gDimension[1]) return ret;
    if (i < 0) {
        i = 0;
    }
    if (k < 0){
        k = 0;
    }
    if (i > gDimension[0]){
        i = gDimension[0] - 1;
    }
    if (k > gDimension[2]){
        k = gDimension[2] - 1;
    }

    int x = i;
    int y = j * (gDimension[0]);
    int z = k * (gDimension[0]) * (gDimension[1]+1);

    return data[x+y+z];
}

vec3 MACGridDataY::worldToLocal(const vec3 &pt) const {
    vec3 ret;
    ret[0] = min(max(0.0f, pt[0]  - CellSize*0.5f), containerBounds[0]);
    ret[1] = min(max(0.0f, pt[1]), containerBounds[1]); //staggering // its already at 0 face of y axis
    ret[2] = min(max(0.0f, pt[2] - CellSize*0.5f), containerBounds[2]);
    return ret;
}

void MACGridDataY::setCell(int &i, int &j, int &k, const float val){
    int x = i;
    int y = j * (gDimension[0]);
    int z = k * (gDimension[0]) * (gDimension[1]+1);

    if(x+y+z < data.size()-1 && x+y+z >= 0)
        data.at(x+y+z) = val;
}

void MACGridDataY::setCellAdd(const int &i, const int &j, const int &k, const float val){
    int x = i;
    int y = j * (gDimension[0]);
    int z = k * (gDimension[0]) * (gDimension[1]+1);

    if(x+y+z < data.size()-1 && x+y+z >= 0)
        data.at(x+y+z) += val;
}

MACGridDataZ::MACGridDataZ()
{
}

void MACGridDataZ::MACGridDataInitialize(){
    containerBounds[0] = CellSize * (gDimension[0]);
    containerBounds[1] = CellSize * gDimension[1];
    containerBounds[2] = CellSize * (gDimension[2] + 1);
    data.resize((gDimension[0]) * gDimension[1] * (gDimension[2] + 1));
    std::fill(data.begin(), data.end(), 0.f);
}

float& MACGridDataZ::operator ()(int i, int j, int k){
    float ret = 0.f;
    if (k < 0 || k > gDimension[2]) return ret;
    if (j < 0) {
        j = 0;
    }
    if (i < 0){
        i = 0;
    }
    if (j > gDimension[1]){
        j = gDimension[1] - 1;
    }
    if (i > gDimension[0]){
        i = gDimension[0] - 1;
    }

    int x = i;
    int y = j * (gDimension[0]);
    int z = k * (gDimension[0]) * (gDimension[1]);

    return data[x+y+z];
}

vec3 MACGridDataZ::worldToLocal(const vec3 &pt) const {
    vec3 ret;
    ret[0] = min(max(0.0f, pt[0] - CellSize*0.5f), containerBounds[0]);
    ret[1] = min(max(0.0f, pt[1] - CellSize*0.5f), containerBounds[1]);
    ret[2] = min(max(0.0f, pt[2]), containerBounds[2]);//staggering // its already at 0 face of z axis
    return ret;
}

void MACGridDataZ::setCell(int &i, int &j, int &k, const float val){
    int x = i;
    int y = j * (gDimension[0]);
    int z = k * (gDimension[0]) * (gDimension[1]);
    data[x+y+z] = val;
}
void MACGridDataZ::setCellAdd(int &i, int &j, int &k, const float val){
    int x = i;
    int y = j * (gDimension[0]);
    int z = k * (gDimension[0]) * (gDimension[1]);
    data[x+y+z] += val;
}
