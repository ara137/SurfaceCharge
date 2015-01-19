/*
Highly Optimized Object-oriented Many-particle Dynamics -- Blue Edition
(HOOMD-blue) Open Source Software License Copyright 2009-2014 The Regents of
the University of Michigan All rights reserved.

HOOMD-blue may contain modifications ("Contributions") provided, and to which
copyright is held, by various Contributors who have granted The Regents of the
University of Michigan the right to modify and/or distribute such Contributions.

You may redistribute, use, and create derivate works of HOOMD-blue, in source
and binary forms, provided you abide by the following conditions:

* Redistributions of source code must retain the above copyright notice, this
list of conditions, and the following disclaimer both in the code and
prominently in any materials provided with the distribution.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions, and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* All publications and presentations based on HOOMD-blue, including any reports
or published results obtained, in whole or in part, with HOOMD-blue, will
acknowledge its use according to the terms posted at the time of submission on:
http://codeblue.umich.edu/hoomd-blue/citations.html

* Any electronic documents citing HOOMD-Blue will link to the HOOMD-Blue website:
http://codeblue.umich.edu/hoomd-blue/

* Apart from the above required attributions, neither the name of the copyright
holder nor the names of HOOMD-blue's contributors may be used to endorse or
promote products derived from this software without specific prior written
permission.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS'' AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND/OR ANY
WARRANTIES THAT THIS SOFTWARE IS FREE OF INFRINGEMENT ARE DISCLAIMED.

IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// **********************
// This code is intended for mimicking the surface charge on aggregated polymers dispersed in water
// First, the polymer center of masses are identified, then a cluster analysis is performed. Then
// the force between the clusters is calculated and set

// inclusion guard
#ifndef _SURFACE_CHARGE_H_
#define _SURFACE_CHARGE_H_

// First, hoomd.h should be included
#include <hoomd/hoomd.h>


class SurfaceCharge : public ForceCompute
    {
    public:
        //! Constructor
        SurfaceCharge(boost::shared_ptr<SystemDefinition> sysdef, boost::shared_ptr<ParticleGroup> group, unsigned int polymer_length);
        //! Function for setting clustering parameters
        void setClusterParams(Scalar);
        //! Function for setting potential parameters
        void setForceParams(Scalar, Scalar, Scalar);
    protected:
        //! Take one timestep forward
        virtual void computeForces(unsigned int timestep);

    private:
        //! Group of particles to apply force to
        boost::shared_ptr<ParticleGroup> m_group;
        //! Internal map of polymers to particle indices
        GPUArray<unsigned int> m_polymer_index_map;
        //! Internal list of polymer center of masses
        GPUArray<Scalar3> m_polymer_com;
        //! Internal list of cluster center of masses
        GPUArray<Scalar3> m_cluster_com;
        //! Internal list of cluster radius of gyration
        GPUArray<Scalar> m_cluster_rg;
        //! Polymer length
        unsigned int m_polymer_length;
        //! Polymer count
        unsigned int m_polymer_count;
        //! Cluster count
        unsigned int m_cluster_count;
        //! Indexes elements in the polymer map
        Index2D m_polymer_indexer;
        //! ClusterIds of polymers
        vector<int> m_cluster_ids;
        //! Neighbor list
        vector< vector<unsigned int> > m_neighbors;
        //! Cluster list storing members of each cluster
        vector< vector<unsigned int> > m_clusters;
        //! Cutoff radius for cluster analysis
        Scalar m_cluster_rcut2;
        //! Parameters for internal potential
        Scalar m_pot_epsilon, m_pot_kappa, m_pot_rcut2;
        //! Parameters for checking whether all parameters have been
        bool m_cluster_params_set, m_force_params_set;

        //! Internal function for remapping the particles to polymers
        void RemapPolymers();
        //! Internal function for calculating the polymer center of masses
        void CalcCenterOfMasses();
        //! Internal function for performing the cluster analysis
        void ClusterAnalysis();
        //! Internal function for calculating the forces between the clusters
        void Forces();
        //! Internal function for recursively assigning cluster ids to the clusters
        void MarkParticleAndNeighbors(unsigned int, unsigned int*);
    };

//! Export the SurfaceCharge class to python
void export_SurfaceCharge();

#endif // _SURFACE_CHARGE_H_

