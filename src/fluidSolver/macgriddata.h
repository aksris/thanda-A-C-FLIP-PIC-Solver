#ifndef MACGRIDDATA_H
#define MACGRIDDATA_H
#include <glm/glm.hpp>
#include <vector>
#include <algorithm>

using namespace glm;

class MACGridData
{
public:
    MACGridData(){};
    MACGridData(ivec3 dimensions, const vec3 &containerBounds, float cellSize);
    MACGridData(const MACGridData &val);

    virtual MACGridData& operator=(const MACGridData& val);
    virtual void MACGridDataInitialize();

    float& operator()(int i, int j, int k);

    virtual void getCell(const vec3& pt, int& i, int& j, int& k);
    int getCellMark(int i, int j, int k);
    int getCellIndex(int i, int j, int k);
    void setCell(int& i, int& j, int& k, const float val);
    void setCellAdd(const int& i,const int& j,const int& k, const float val);
    void setCellMark(int &i, int &j, int &k, const int val, bool mark);
    virtual vec3 worldToLocal(const vec3& pt) const;
    virtual float interpolate(const vec3& pt);
    vec3 containerBounds;
    std::vector<float> data;
    std::vector<int> mData;
    float CellSize;
    ivec3 resolution;

protected:
};

class MACGridDataX : public MACGridData
{
public:
    MACGridDataX():MACGridData(){};
    MACGridDataX(ivec3 dimensions, const vec3& containerBounds, float cellSize) : MACGridData(dimensions, containerBounds, cellSize){};
    virtual void MACGridDataInitialize();
    virtual vec3 worldToLocal(const vec3& pt) const;
};

class MACGridDataY : public MACGridData
{
public:
    MACGridDataY():MACGridData(){};
    MACGridDataY(ivec3 dimensions, const vec3& containerBounds, float cellSize) : MACGridData(dimensions, containerBounds, cellSize){};
    virtual void MACGridDataInitialize();
    virtual vec3 worldToLocal(const vec3& pt) const;
};

class MACGridDataZ : public MACGridData
{
public:
    MACGridDataZ():MACGridData(){};
    MACGridDataZ(ivec3 dimensions, const vec3& containerBounds, float cellSize) : MACGridData(dimensions, containerBounds, cellSize){};
    virtual void MACGridDataInitialize();
    virtual vec3 worldToLocal(const vec3& pt) const;
};


#endif // MACGRIDDATA_H
