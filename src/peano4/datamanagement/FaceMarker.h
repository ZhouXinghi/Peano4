// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once


#include "peano4/utils/Globals.h"
#include "peano4/grid/GridTraversalEvent.h"
#include "tarch/la/Vector.h"


#include <bitset>


namespace peano4 {
  namespace datamanagement {
    struct FaceMarker;
  }
}


std::ostream& operator<<( std::ostream& out, const peano4::datamanagement::FaceMarker& marker );


/**
 * Provide information about selected face
 *
 * This routine is passed into any observer/action set routine which
 * refers to a face, so routines like touchFaceFirstTime(). It provides
 * you with information about the particular face as well as their location
 * within the mesh.
 *
 * @see VertexMarker for an analogous description of the vertex handling.
 * @see FaceEnumerator for a discussion of the ordering of the faces.
 */
struct peano4::datamanagement::FaceMarker {
  private:
    /**
     * Centre of the underlying cell
     */
    tarch::la::Vector<Dimensions,double>  _cellCentre;

    /**
     * Size of the underlying cell
     */
    tarch::la::Vector<Dimensions,double>  _h;

    std::bitset<2*Dimensions>             _hasBeenRefined;
    std::bitset<2*Dimensions>             _willBeRefined;
    std::bitset<2*Dimensions>             _isLocal;

    bool                                  _cellIsLocal;

    int _select;


    /**
     * Entries from (0,1,2). (0,0) or (0,0,0) is the left, bottom cell.
     */
    tarch::la::Vector<Dimensions,int>  _relativePositionOfCellWithinFatherCell;
  public:
    /**
     * The derivation of _isLocal and _isRefined is very similar to
     * peano4::grid::Spacetree::getFaceType().
     */
    FaceMarker(const peano4::grid::GridTraversalEvent& event, int select=0);

    /**
     * Selects a face within a cell, i.e. now the marker knows to which
     * face it corresponds. After that, the routine returns a this
     * reference.
     *
     * @param face Number from 0 to 2d-1
     */
    FaceMarker& select(int face);

    /**
     * Return the number of the face.
     *
     * If you use the modulo Dimensions, then you implicitly know along which
     * coordinate axis the corresponding face normal points.
     *
     * @see FaceEnumerator for a discussion of the face ordering.
     *
     * @return Number from 0 to 2d-1
     */
    int getSelectedFaceNumber() const;

    /**
     * Center of a particular face with respective reference cell.
     */
    tarch::la::Vector<Dimensions,double> x(int i) const;

    /**
     * A marker is always tied to one particular face. You can still
     * get the data of all faces of a cell (see x(int)), but this
     * particular function gives you the x coordinate with the centre
     * of the currently selected face.
     *
     * @see select(int)
     */
    tarch::la::Vector<Dimensions,double> x() const;

    /**
     * This operation gives you the normal of a certain face.
     *
     * The orientation of the normal depends on the context in which
     * the face marker is used: Every operation within Peano is
     * always triggered from a cell point of view. The normal corresponds
     * to this view. If you run through a mesh twice, and if you twice
     * hit the same face, the normal thus might be mirrored as a
     * face marker might be constructed from the other adjacent cell.
     *
     * We do enumerate the faces as follows: The first face is the one
     * whose normal is parallel to first coordinate axis. The second
     * face is the one whose normal is parallel to the second axis. The
     * d+1th face is the one parallel to the very first face. The very
     * first face however is closer to the origin in the unit cube.
     *
     * @image html FaceEnumerator.png
     */
    tarch::la::Vector<Dimensions,double> normal(int i) const;
    tarch::la::Vector<Dimensions,double> normal() const;

    /**
     * This operation gives you the outer normal of a cell. Different to
     * normal(int), the normal however is oriented along the domain
     * boundaries, i.e. it always points outside from the local domain.
     */
    tarch::la::Vector<Dimensions,double> outerNormal(int i) const;
    tarch::la::Vector<Dimensions,double> outerNormal() const;

    /**
     * Size of the underlying cell
     */
    tarch::la::Vector<Dimensions,double>  h() const;

    std::string toString() const;

    bool hasBeenRefined() const;
    bool hasBeenRefined(int i) const;

    bool willBeRefined() const;
    bool willBeRefined(int i) const;

    bool isLocal() const;
    bool isLocal(int i) const;

    /**
     * Return relative position within father cell. The result is from
     * (0,1,2,3) x (0,1,2) x (0,1,2) or (0,1,2) x (0,1,2,3) x (0,1,2) or
     * (0,1,2) x (0,1,2) x (0,1,2,3).
     */
    tarch::la::Vector<Dimensions,int>  getRelativePositionWithinFatherCell() const;
    tarch::la::Vector<Dimensions,int>  getRelativePositionWithinFatherCell(int i) const;

    /**
     * Return relative position of a subface within its father face. The
     * result is from (0,1,2)^(Dimensions-1)
     */
    tarch::la::Vector<Dimensions-1,int>  getRelativePositionWithinFatherFace() const;
    
    /**
     * The term patch here refers to a 3x3 or 3x3x3 refinement. The routine
     * determines whether the face is on the boundary of this patch of not.
     */
    bool isInteriorFaceWithinPatch() const;
};


