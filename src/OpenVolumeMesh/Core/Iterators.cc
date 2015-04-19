/*===========================================================================*\
 *                                                                           *
 *                            OpenVolumeMesh                                 *
 *        Copyright (C) 2011 by Computer Graphics Group, RWTH Aachen         *
 *                        www.openvolumemesh.org                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of OpenVolumeMesh.                                     *
 *                                                                           *
 *  OpenVolumeMesh is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenVolumeMesh is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenVolumeMesh.  If not,                              *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision: 267 $                                                         *
 *   $Date: 2014-05-02 12:21:27 +0200 (Fri, 02 May 2014) $                    *
 *   $LastChangedBy: kremer $                                                *
 *                                                                           *
\*===========================================================================*/

#include <algorithm>
#ifndef NDEBUG
#include <iostream>
#endif
#include <set>

#include "Iterators.hh"
#include "TopologyKernel.hh"
#include "OpenVolumeMeshHandle.hh"

namespace OpenVolumeMesh {

//================================================================================================
// VertexOHalfEdgeIter
//================================================================================================


VertexOHalfEdgeIter::VertexOHalfEdgeIter(const VertexHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps),
cur_index_(0) {

	if(!_mesh->has_vertex_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

	if((unsigned int)_ref_h.idx() >= BaseIter::mesh()->outgoing_hes_per_vertex_.size()) {
		BaseIter::valid(false);
	}

	if(BaseIter::valid()) {
		if((unsigned int)cur_index_ >= BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h.idx()].size()) {
			BaseIter::valid(false);
		}
	}

	if(BaseIter::valid()) {
		BaseIter::cur_handle((
				BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h.idx()])[cur_index_]);
	}
}


VertexOHalfEdgeIter& VertexOHalfEdgeIter::operator--() {

    size_t n_outgoing_halfedges = BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle().idx()].size();

    if (cur_index_ == 0) {
        cur_index_ = n_outgoing_halfedges-1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle((BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle().idx()])[cur_index_]);

	return *this;
}


VertexOHalfEdgeIter& VertexOHalfEdgeIter::operator++() {

    size_t n_outgoing_halfedges = BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle().idx()].size();

    ++cur_index_;

    if (cur_index_ == n_outgoing_halfedges) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle((BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle().idx()])[cur_index_]);

	return *this;
}

////================================================================================================
//// HalfEdgeHalfFaceIter
////================================================================================================


HalfEdgeHalfFaceIter::HalfEdgeHalfFaceIter(const HalfEdgeHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps),
cur_index_(0) {

	if(!_mesh->has_edge_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

	if((unsigned int)_ref_h.idx() >= BaseIter::mesh()->incident_hfs_per_he_.size()) {
		BaseIter::valid(false);
	}

	if(BaseIter::valid()) {
		if((unsigned int)cur_index_ >= BaseIter::mesh()->incident_hfs_per_he_[_ref_h.idx()].size()) {
			BaseIter::valid(false);
		}
	}

	if(BaseIter::valid()) {
		BaseIter::cur_handle((
				BaseIter::mesh()->incident_hfs_per_he_[_ref_h.idx()])[cur_index_]);
	}
}


HalfEdgeHalfFaceIter& HalfEdgeHalfFaceIter::operator--() {

    size_t n_outgoing_halffaces = BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle().idx()].size();

    if (cur_index_ == 0) {
        cur_index_ = n_outgoing_halffaces-1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle((BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle().idx()])[cur_index_]);

    return *this;
}


HalfEdgeHalfFaceIter& HalfEdgeHalfFaceIter::operator++() {


    size_t n_outgoing_halffaces = BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle().idx()].size();

    ++cur_index_;

    if (cur_index_ == n_outgoing_halffaces) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle((BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle().idx()])[cur_index_]);

    return *this;
}

////================================================================================================
//// VertexCellIter
////================================================================================================

VertexCellIter::VertexCellIter(const VertexHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

	if(!_mesh->has_full_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

    if((unsigned int)_ref_h.idx() >= BaseIter::mesh()->outgoing_hes_per_vertex_.size()) {
        BaseIter::valid(false);
        return;
    }

    // Build up cell list
    const std::vector<HalfEdgeHandle>& incidentHalfedges = BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h.idx()];
    for(std::vector<HalfEdgeHandle>::const_iterator it = incidentHalfedges.begin(); it != incidentHalfedges.end(); ++it) {

    	if(*it < 0 || (unsigned int)it->idx() >= BaseIter::mesh()->incident_hfs_per_he_.size()) continue;
        const std::vector<HalfFaceHandle>& incidentHalfFaces = BaseIter::mesh()->incident_hfs_per_he_[it->idx()];

    	for(std::vector<HalfFaceHandle>::const_iterator hf_it = incidentHalfFaces.begin();
                hf_it != incidentHalfFaces.end(); ++hf_it) {
    		if((unsigned int)hf_it->idx() < BaseIter::mesh()->incident_cell_per_hf_.size()) {
    			CellHandle c_idx = BaseIter::mesh()->incident_cell_per_hf_[hf_it->idx()];
    			if(c_idx != TopologyKernel::InvalidCellHandle)
                    cells_.push_back(c_idx);
    		}
    	}
    }
    // Remove all duplicate entries
    std::sort(cells_.begin(), cells_.end());
    cells_.resize(std::unique(cells_.begin(), cells_.end()) - cells_.begin());

    cur_index_ = 0;
    BaseIter::valid(cells_.size()>0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(cells_[cur_index_]);
    }
}

VertexCellIter& VertexCellIter::operator--() {

    if(cur_index_ == 0) {
        cur_index_ = cells_.size()-1;
        --lap_;
        if (lap_ < 0) {
            BaseIter::valid(false);
        }
    } else {
        --cur_index_;
    }

    BaseIter::cur_handle(cells_[cur_index_]);
    return *this;
}


VertexCellIter& VertexCellIter::operator++() {

    ++cur_index_;
    if(cur_index_ == cells_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_) {
            BaseIter::valid(false);
        }
	}
    BaseIter::cur_handle(cells_[cur_index_]);
	return *this;
}

////================================================================================================
//// HalfEdgeCellIter
////================================================================================================


HalfEdgeCellIter::HalfEdgeCellIter(const HalfEdgeHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps),
cur_index_(0) {

	if(!_mesh->has_edge_bottom_up_incidences() || !_mesh->has_face_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

    if((unsigned int)_ref_h.idx() >= BaseIter::mesh()->incident_hfs_per_he_.size()) {

        BaseIter::valid(false);
        return;
    }
    if((unsigned int)cur_index_ >= BaseIter::mesh()->incident_hfs_per_he_[_ref_h.idx()].size()) {

    	BaseIter::valid(false);
    	return;
    }
    if((unsigned int)((BaseIter::mesh()->incident_hfs_per_he_[_ref_h.idx()])[cur_index_]).idx() >=
    		BaseIter::mesh()->incident_cell_per_hf_.size()) {

    	BaseIter::valid(false);
    	return;
    }

    // collect cell handles
    const std::vector<HalfFaceHandle>& incidentHalffaces = BaseIter::mesh()->incident_hfs_per_he_[_ref_h.idx()];
    std::set<CellHandle> cells;
    for (unsigned int i = 0; i < incidentHalffaces.size(); ++i)
    {
        CellHandle ch = getCellHandle(i);
        if (ch.is_valid()) {
            if(cells.count(ch) == 0) {
                cells_.push_back(ch);
            }
            cells.insert(ch);
        }
    }

    BaseIter::valid(cells_.size() > 0);

    if(BaseIter::valid()) {
        BaseIter::cur_handle(cells_[cur_index_]);
    }
}


HalfEdgeCellIter& HalfEdgeCellIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_ = cells_.size()-1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle(cells_[cur_index_]);
    return *this;
}


HalfEdgeCellIter& HalfEdgeCellIter::operator++() {

    ++cur_index_;

    if (cur_index_ == cells_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle(cells_[cur_index_]);
    return *this;
}

CellHandle HalfEdgeCellIter::getCellHandle(int _cur_index) const
{
    const std::vector<HalfFaceHandle>& halffacehandles = BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle().idx()];
    HalfFaceHandle currentHalfface = halffacehandles[_cur_index];
    if(!currentHalfface.is_valid()) return CellHandle(-1);
    CellHandle cellhandle = BaseIter::mesh()->incident_cell_per_hf_[currentHalfface.idx()];
    return cellhandle;
}


////================================================================================================
//// CellVertexIter
////================================================================================================


CellVertexIter::CellVertexIter(const CellHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    std::vector<HalfFaceHandle>::const_iterator hf_iter = BaseIter::mesh()->cell(_ref_h).halffaces().begin();
    for(; hf_iter != BaseIter::mesh()->cell(_ref_h).halffaces().end(); ++hf_iter) {
        const OpenVolumeMeshFace& halfface = BaseIter::mesh()->halfface(*hf_iter);
        const std::vector<HalfEdgeHandle>& hes = halfface.halfedges();
        for(std::vector<HalfEdgeHandle>::const_iterator he_iter = hes.begin(); he_iter != hes.end(); ++he_iter) {
            incident_vertices_.push_back(BaseIter::mesh()->halfedge(*he_iter).to_vertex());
        }
    }

    // Remove all duplicate entries
    std::sort(incident_vertices_.begin(), incident_vertices_.end());
    incident_vertices_.resize(std::unique(incident_vertices_.begin(), incident_vertices_.end()) - incident_vertices_.begin());

    cur_index_ = 0;
    BaseIter::valid(incident_vertices_.size() > 0);

    if(BaseIter::valid()) {
        BaseIter::cur_handle(incident_vertices_[cur_index_]);
    }
}


CellVertexIter& CellVertexIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_ = incident_vertices_.size()-1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle(incident_vertices_[cur_index_]);
    return *this;
}


CellVertexIter& CellVertexIter::operator++() {

    ++cur_index_;
    if (cur_index_ == incident_vertices_.size()){
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle(incident_vertices_[cur_index_]);
	return *this;
}

////================================================================================================
//// CellCellIter
////================================================================================================


CellCellIter::CellCellIter(const CellHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    if(!_mesh->has_face_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

	std::vector<HalfFaceHandle>::const_iterator hf_iter = BaseIter::mesh()->cell(_ref_h).halffaces().begin();
    std::vector<HalfFaceHandle>::const_iterator hf_end  = BaseIter::mesh()->cell(_ref_h).halffaces().end();
    for(; hf_iter != hf_end; ++hf_iter) {

		HalfFaceHandle opp_hf = BaseIter::mesh()->opposite_halfface_handle(*hf_iter);
		CellHandle ch = BaseIter::mesh()->incident_cell_per_hf_[opp_hf.idx()];
		if(ch != TopologyKernel::InvalidCellHandle) {
            adjacent_cells_.push_back(ch);
		}
	}

    // Remove all duplicate entries
    std::sort(adjacent_cells_.begin(), adjacent_cells_.end());
    adjacent_cells_.resize(std::unique(adjacent_cells_.begin(), adjacent_cells_.end()) - adjacent_cells_.begin());

    cur_index_ = 0;
    BaseIter::valid(adjacent_cells_.size()>0);
	if(BaseIter::valid()) {
        BaseIter::cur_handle(adjacent_cells_[cur_index_]);
	}
}


CellCellIter& CellCellIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_  = adjacent_cells_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle(adjacent_cells_[cur_index_]);
    return *this;
}


CellCellIter& CellCellIter::operator++() {

    ++cur_index_;
    if (cur_index_ == adjacent_cells_.size())
    {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(adjacent_cells_[cur_index_]);
    return *this;
}

////================================================================================================
//// HalfFaceVertexIter
////================================================================================================


HalfFaceVertexIter::HalfFaceVertexIter(const HalfFaceHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    if(!_ref_h.is_valid()) return;

    const OpenVolumeMeshFace& halfface = _mesh->halfface(_ref_h);
    const std::vector<HalfEdgeHandle>& hes = halfface.halfedges();
    for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin();
            he_it != hes.end(); ++he_it) {
        vertices_.push_back(_mesh->halfedge(*he_it).from_vertex());
    }

    cur_index_ = 0;

    BaseIter::valid(vertices_.size() > 0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(vertices_[cur_index_]);
    }
}


HalfFaceVertexIter& HalfFaceVertexIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_  = vertices_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle(vertices_[cur_index_]);
    return *this;
}


HalfFaceVertexIter& HalfFaceVertexIter::operator++() {

    ++cur_index_;
    if (cur_index_ == vertices_.size())
    {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(vertices_[cur_index_]);
    return *this;
}

//================================================================================================
// BoundaryHalfFaceHalfFaceIter
//================================================================================================

BoundaryHalfFaceHalfFaceIter::BoundaryHalfFaceHalfFaceIter(const HalfFaceHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    if(!_mesh->has_face_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

    // Go over all incident halfedges
//    const std::vector<HalfEdgeHandle> halfedges = _mesh->halfface(_ref_h).halfedges();
    const OpenVolumeMeshFace& halfface = _mesh->halfface(_ref_h);
    const std::vector<HalfEdgeHandle>& halfedges = halfface.halfedges();
    for(std::vector<HalfEdgeHandle>::const_iterator he_it = halfedges.begin();
            he_it != halfedges.end(); ++he_it) {

        // Get outside halffaces
        OpenVolumeMesh::HalfEdgeHalfFaceIter hehf_it = _mesh->hehf_iter(_mesh->opposite_halfedge_handle(*he_it));
        for(; hehf_it.valid(); ++hehf_it) {

            if(_mesh->is_boundary(*hehf_it)) {
                neighbor_halffaces_.push_back(*hehf_it);
                common_edges_.push_back(_mesh->edge_handle(*he_it));
            }
        }
    }

    cur_index_ = 0;
    BaseIter::valid(neighbor_halffaces_.size() > 0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(neighbor_halffaces_[cur_index_]);
    }
}


BoundaryHalfFaceHalfFaceIter& BoundaryHalfFaceHalfFaceIter::operator--() {

    if (cur_index_ == 0)
    {
        cur_index_ = neighbor_halffaces_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else{
        --cur_index_;
    }

    BaseIter::cur_handle(neighbor_halffaces_[cur_index_]);
    return *this;
}


BoundaryHalfFaceHalfFaceIter& BoundaryHalfFaceHalfFaceIter::operator++() {

    ++cur_index_;
    if (cur_index_ == neighbor_halffaces_.size())
    {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle(neighbor_halffaces_[cur_index_]);
    return *this;
}

////================================================================================================
//// BoundaryFaceIter
////================================================================================================


BoundaryFaceIter::BoundaryFaceIter(const TopologyKernel* _mesh) :
BaseIter(_mesh),
bf_it_(_mesh->faces_begin()) {

	if(!_mesh->has_face_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

	while(bf_it_ != BaseIter::mesh()->faces_end() &&
	        !BaseIter::mesh()->is_boundary(*bf_it_)) {
	    ++bf_it_;
	}
	BaseIter::valid(bf_it_ != BaseIter::mesh()->faces_end());
	if(BaseIter::valid()) {
		BaseIter::cur_handle(*bf_it_);
	}
}


BoundaryFaceIter& BoundaryFaceIter::operator--() {

    --bf_it_;
    while(bf_it_ >= BaseIter::mesh()->faces_begin() &&
            !BaseIter::mesh()->is_boundary(*bf_it_)) {
        --bf_it_;
    }
	if(bf_it_ >= BaseIter::mesh()->faces_begin()) {
		BaseIter::cur_handle(*bf_it_);
	} else {
		BaseIter::valid(false);
	}
	return *this;
}


BoundaryFaceIter& BoundaryFaceIter::operator++() {

	++bf_it_;
	while(bf_it_ != BaseIter::mesh()->faces_end() &&
            !BaseIter::mesh()->is_boundary(*bf_it_)) {
        ++bf_it_;
    }
	if(bf_it_ != BaseIter::mesh()->faces_end()) {
		BaseIter::cur_handle(*bf_it_);
	} else {
		BaseIter::valid(false);
	}
	return *this;
}

////================================================================================================
//// VertexIter
////================================================================================================


VertexIter::VertexIter(const TopologyKernel* _mesh, const VertexHandle& _vh) :
BaseIter(_mesh, _vh),
cur_index_(_vh.idx()) {

	if((unsigned int)cur_index_ >= BaseIter::mesh()->n_vertices()) {
		BaseIter::valid(false);
	}
	if(BaseIter::valid()) {
		BaseIter::cur_handle(VertexHandle(cur_index_));
	}
}


VertexIter& VertexIter::operator--() {

    --cur_index_;
    if(cur_index_ < 0) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(VertexHandle(cur_index_));
    return *this;
}


VertexIter& VertexIter::operator++() {

	++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->n_vertices()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(VertexHandle(cur_index_));
	return *this;
}

////================================================================================================
//// EdgeIter
////================================================================================================


EdgeIter::EdgeIter(const TopologyKernel* _mesh, const EdgeHandle& _eh) :
BaseIter(_mesh, _eh),
cur_index_(_eh.idx()) {

	if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size()) {
		BaseIter::valid(false);
	}
	if(BaseIter::valid()) {
		BaseIter::cur_handle(EdgeHandle(cur_index_));
	}
}


EdgeIter& EdgeIter::operator--() {

	--cur_index_;
	if(cur_index_ < 0) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(EdgeHandle(cur_index_));
	return *this;
}


EdgeIter& EdgeIter::operator++() {

	++cur_index_;
	if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size()) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(EdgeHandle(cur_index_));
	return *this;
}

////================================================================================================
//// HalfEdgeIter
////================================================================================================


HalfEdgeIter::HalfEdgeIter(const TopologyKernel* _mesh, const HalfEdgeHandle& _heh) :
BaseIter(_mesh, _heh),
cur_index_(_heh.idx()) {

	if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size() * 2) {
		BaseIter::valid(false);
	}
	if(BaseIter::valid()) {
		BaseIter::cur_handle(HalfEdgeHandle(cur_index_));
	}
}


HalfEdgeIter& HalfEdgeIter::operator--() {

	--cur_index_;
	if(cur_index_ < 0) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(HalfEdgeHandle(cur_index_));
	return *this;
}


HalfEdgeIter& HalfEdgeIter::operator++() {

	++cur_index_;
	if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size() * 2) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(HalfEdgeHandle(cur_index_));
	return *this;
}

////================================================================================================
//// FaceIter
////================================================================================================


FaceIter::FaceIter(const TopologyKernel* _mesh, const FaceHandle& _fh) :
BaseIter(_mesh, _fh),
cur_index_(_fh.idx()) {

	if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size()) {
		BaseIter::valid(false);
	}
	if(BaseIter::valid()) {
		BaseIter::cur_handle(FaceHandle(cur_index_));
	}
}


FaceIter& FaceIter::operator--() {

	--cur_index_;
	if(cur_index_ < 0) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(FaceHandle(cur_index_));
	return *this;
}


FaceIter& FaceIter::operator++() {

	++cur_index_;
	if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size()) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(FaceHandle(cur_index_));
	return *this;
}

////================================================================================================
//// HalfFaceIter
////================================================================================================


HalfFaceIter::HalfFaceIter(const TopologyKernel* _mesh, const HalfFaceHandle& _hfh) :
BaseIter(_mesh, _hfh),
cur_index_(_hfh.idx()) {

	if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size() * 2) {
		BaseIter::valid(false);
	}
	if(BaseIter::valid()) {
		BaseIter::cur_handle(HalfFaceHandle(cur_index_));
	}
}


HalfFaceIter& HalfFaceIter::operator--() {

	--cur_index_;
	if(cur_index_ < 0) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(HalfFaceHandle(cur_index_));
	return *this;
}


HalfFaceIter& HalfFaceIter::operator++() {

	++cur_index_;
	if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size() * 2) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(HalfFaceHandle(cur_index_));
	return *this;
}

////================================================================================================
//// CellIter
////================================================================================================


CellIter::CellIter(const TopologyKernel* _mesh, const CellHandle& _ch) :
BaseIter(_mesh, _ch),
cur_index_(_ch.idx()) {

	if((unsigned int)cur_index_ >= BaseIter::mesh()->cells_.size()) {
		BaseIter::valid(false);
	}
	if(BaseIter::valid()) {
		BaseIter::cur_handle(CellHandle(cur_index_));
	}
}


CellIter& CellIter::operator--() {

	--cur_index_;
	if(cur_index_ < 0) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(CellHandle(cur_index_));
	return *this;
}


CellIter& CellIter::operator++() {

	++cur_index_;
	if((unsigned int)cur_index_ >= BaseIter::mesh()->cells_.size()) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(CellHandle(cur_index_));
	return *this;
}

} // Namespace OpenVolumeMesh
