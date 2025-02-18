///////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2012-2016 DreamWorks Animation LLC
//
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
//
// Redistributions of source code must retain the above copyright
// and license notice and the following restrictions and disclaimer.
//
// *     Neither the name of DreamWorks Animation nor the names of
// its contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// IN NO EVENT SHALL THE COPYRIGHT HOLDERS' AND CONTRIBUTORS' AGGREGATE
// LIABILITY FOR ALL CLAIMS REGARDLESS OF THEIR BASIS EXCEED US$250.00.
//
///////////////////////////////////////////////////////////////////////////

#include <cppunit/extensions/HelperMacros.h>
#include <openvdb/tools/MaskToLevelSet.h>

#define ASSERT_DOUBLES_EXACTLY_EQUAL(expected, actual) \
    CPPUNIT_ASSERT_DOUBLES_EQUAL((expected), (actual), /*tolerance=*/0.0);

class MaskToLevelSet: public CppUnit::TestCase
{
public:
    CPPUNIT_TEST_SUITE(MaskToLevelSet);
    CPPUNIT_TEST(testConversion);
    CPPUNIT_TEST_SUITE_END();

    void testConversion();
};

CPPUNIT_TEST_SUITE_REGISTRATION(MaskToLevelSet);

void MaskToLevelSet::testConversion()
{
    using namespace openvdb;
    typedef tree::Tree4<ValueMask, 5, 4, 3>::Type MyMaskTree;
    typedef Grid<MyMaskTree>                      MyMaskGrid;

    typedef tree::Tree4<float, 5, 4, 3>::Type     MyFloatTree;
    typedef Grid<MyFloatTree>                     MyFloatGrid;

    const float voxelSize = 0.1f;

    MyMaskGrid maskGrid( false );
    maskGrid.setTransform( math::Transform::createLinearTransform( voxelSize ) );

    // Define active region
    maskGrid.fill( CoordBBox( Coord(0), Coord(7)), true );
    maskGrid.tree().voxelizeActiveTiles();
    
    MyFloatGrid::Ptr sdfGrid = tools::maskToLevelSet( maskGrid );

    CPPUNIT_ASSERT( sdfGrid.get() != NULL );
    CPPUNIT_ASSERT( !sdfGrid->empty() );
    CPPUNIT_ASSERT_EQUAL( int( GRID_LEVEL_SET ), int( sdfGrid->getGridClass() ) );
    ASSERT_DOUBLES_EXACTLY_EQUAL( voxelSize, sdfGrid->voxelSize()[0] );

    // test inside coord value
    CPPUNIT_ASSERT( sdfGrid->tree().getValue( Coord(3,3,3)) < 0.0f );

    // test outside coord value
    CPPUNIT_ASSERT( sdfGrid->tree().getValue( Coord(10,10,10)) > 0.0f );
}


// Copyright (c) 2012-2016 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
