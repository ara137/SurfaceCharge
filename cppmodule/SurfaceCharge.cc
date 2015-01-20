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

// we need to include boost.python in order to export SurfaceCharge to python
#include <boost/python.hpp>
using namespace boost::python;

// we need to include boost.bind for GPUWorker execution
#include <boost/bind.hpp>
using namespace boost;

#include "SurfaceCharge.h"

// ********************************
// here follows the code for SurfaceCharge on the CPU

/*! \param sysdef System to zero the velocities of
*/
SurfaceCharge::SurfaceCharge(boost::shared_ptr<SystemDefinition> sysdef, boost::shared_ptr<ParticleGroup> group, unsigned int polymer_length)
        : ForceCompute(sysdef)
    {
    // Assign group and resize polymer map
    m_group = group;
    m_polymer_index_map.resize(group->getNumMembers());
    m_polymer_length = polymer_length;

    // Set flag that parameters need to be set
    m_cluster_params_set = false;
    m_force_params_set = false;
    
    // Check for stray monomers
    if (m_group->getNumMembers()%m_polymer_length != 0)
        {
        std::cerr<<"SurfaceCharge: Number of particles is not an integer multiple of monomers/chain!"<<std::endl;
        }

    // Set polymer count
    m_polymer_count = m_group->getNumMembers()/m_polymer_length;
    // Resize center of masses array
    m_polymer_com.resize(m_polymer_count);
    // Resize clusterIds array
    m_cluster_ids.resize(m_polymer_count);

    // Initialize indexer, where each row consists of one polymer
    m_polymer_indexer = Index2D(m_polymer_length, m_polymer_count);
    // Map has to be created when the object is constructed
    RemapPolymers();
    }

/*! Map the particle indices to polymers
    For the time being, use a fixed length to separate polymers from each other
*/
void SurfaceCharge::RemapPolymers()
    {
    assert(m_pdata);
    assert(m_polymer_index_map);
    ArrayHandle<unsigned int> h_polymer_index_map(m_polymer_index_map, access_location::host, access_mode::readwrite);
    ArrayHandle<unsigned int> h_tag(m_pdata->getTags(), access_location::host, access_mode::read);

    //std::cout<<"SurfaceCharge: Remapping Polymers!"<<std::endl;
    // Update the mapping
    for (unsigned int i=0; i<m_group->getNumMembers(); i++)
        {
        // get the index for the current group member
        const unsigned int idx = m_group->getMemberIndex(i);
        const unsigned int polymer_id = h_tag.data[idx]/m_polymer_length;
        const unsigned int monomer_id = h_tag.data[idx]%m_polymer_length;
        h_polymer_index_map.data[m_polymer_indexer(monomer_id, polymer_id)] = idx;
        }
    }

/*! Calculate the polymer center of masses using the established map
*/
void SurfaceCharge::CalcCenterOfMasses()
    {
    assert(m_pdata);
    assert(m_polymer_index_map);
    assert(m_polymer_com);
    
    // Acquire particle data and simulation box
    ArrayHandle<unsigned int> h_polymer_index_map(m_polymer_index_map, access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_polymer_com(m_polymer_com, access_location::host, access_mode::readwrite);
    const BoxDim& box = m_pdata->getBox();    

    for (unsigned int i=0; i<m_polymer_count; i++)
        {
        Scalar3 com = make_scalar3(0.0, 0.0, 0.0);

        // index of the first bead
        const unsigned int idx_first = h_polymer_index_map.data[m_polymer_indexer(0, i)]; 

        for (unsigned int j=1; j<m_polymer_length; j++)
            {
            const unsigned int idx = h_polymer_index_map.data[m_polymer_indexer(j, i)];
            
            Scalar3 dvec = make_scalar3(h_pos.data[idx].x - h_pos.data[idx_first].x,
                                        h_pos.data[idx].y - h_pos.data[idx_first].y,
                                        h_pos.data[idx].z - h_pos.data[idx_first].z);
            dvec = box.minImage(dvec);
            com.x += dvec.x;
            com.y += dvec.y;
            com.z += dvec.z;
            }
            
        com.x /= m_polymer_length;
        com.y /= m_polymer_length;
        com.z /= m_polymer_length;
            
        com.x += h_pos.data[idx_first].x;
        com.y += h_pos.data[idx_first].y;
        com.z += h_pos.data[idx_first].z;

        int3 img_dummy;
        box.wrap(com, img_dummy);
        h_polymer_com.data[i] = com;
        //std::cout<<"SurfaceCharge: polymer "<<i<<", center of mass: "<<com.x<<", "<<com.y<<", "<<com.z<<std::endl;
        }
    }

/*! Perform Cluster Analysis on the identified polymer center of masses
*/
void SurfaceCharge::ClusterAnalysis()
    {
    assert(m_pdata);
    assert(m_polymer_com);
    assert(m_polymer_index_map);
    assert(m_cluster_com);
    assert(m_cluster_rg);

    // Acquire handle to required data and simulation box
    ArrayHandle<Scalar3> h_polymer_com(m_polymer_com, access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_polymer_index_map(m_polymer_index_map, access_location::host, access_mode::read);
    const BoxDim& box = m_pdata->getBox();

    // Reset the neighborlist and reserve sufficient space in the list
    m_neighbors.assign(m_polymer_count, vector<unsigned int>(0));

    for (unsigned int i=0; i<m_polymer_count; i++)
        {
        m_neighbors[i].reserve(m_polymer_count);
        }

    // Identify neighbors and add them to the lists
    for (unsigned int i=0; i<m_polymer_count; i++)
        {
        for (unsigned int j=i+1; j<m_polymer_count; j++)
            {
            Scalar3 dvec = make_scalar3(h_polymer_com.data[j].x - h_polymer_com.data[i].x,
                                        h_polymer_com.data[j].y - h_polymer_com.data[i].y,
                                        h_polymer_com.data[j].z - h_polymer_com.data[i].z);
            dvec = box.minImage(dvec);
            const Scalar dr2 = dot(dvec, dvec);

            if (dr2 < m_cluster_rcut2)
                {
                m_neighbors[i].push_back(j);
                m_neighbors[j].push_back(i);
                }
            }   
        }

    // Group polymers to clusters and assign unique cluster ids
    unsigned int current_cluster_id = 0;
    // Reset old cluster ids (-1 refers to unassigned polymer)
    for (unsigned int i=0; i<m_polymer_count; i++)
        {
        m_cluster_ids[i] = -1;
        }
    
    // Now go through all polymers, and recursively assign all neighbors a cluster id, until no more neighbors are found
    for (unsigned int i=0; i<m_polymer_count; i++)
        {
            if ((m_cluster_ids[i] == -1) && (i > 0))
                {
                current_cluster_id++;
                }
            MarkParticleAndNeighbors(i, &current_cluster_id);
        }
    
    // Resize internal cluster arrays
    m_cluster_count = current_cluster_id+1;
    m_clusters.assign(m_cluster_count, vector<unsigned int>(0));
    m_cluster_com.resize(m_cluster_count);
    m_cluster_rg.resize(m_cluster_count);

    // Update cluster list
    for (unsigned int i=0; i<m_polymer_count; i++)
        {
        m_clusters[m_cluster_ids[i]].push_back(i);
        }

    // Calculate cluster center of masses of and its radius of gyration
    // Acquire handle to respetive lists
    ArrayHandle<Scalar3> h_cluster_com(m_cluster_com, access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar> h_cluster_rg(m_cluster_rg, access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);

    // Compute center of mass
    for (unsigned int i=0; i<m_cluster_count; i++)
        {
        Scalar3 com = make_scalar3(0.0, 0.0, 0.0);

        // index of the first polymer in the cluster
        const unsigned int idx_first = m_clusters[i][0];

        for (unsigned int j=1; j<m_clusters[i].size(); j++)
            {
            const unsigned int idx = m_clusters[i][j];

            Scalar3 dvec = make_scalar3(h_polymer_com.data[idx].x - h_polymer_com.data[idx_first].x,
                                        h_polymer_com.data[idx].y - h_polymer_com.data[idx_first].y,
                                        h_polymer_com.data[idx].z - h_polymer_com.data[idx_first].z);
            dvec = box.minImage(dvec);
            com.x += dvec.x;
            com.y += dvec.y;
            com.z += dvec.z;
            }
        
        com.x /= m_clusters[i].size();
        com.y /= m_clusters[i].size();
        com.z /= m_clusters[i].size();

        com.x += h_polymer_com.data[idx_first].x;
        com.y += h_polymer_com.data[idx_first].y;
        com.z += h_polymer_com.data[idx_first].z;

        int3 img_dummy;
        box.wrap(com, img_dummy);
        h_cluster_com.data[i] = com;
        }

    // Computer radius of gyration
    for (unsigned int i=0; i<m_cluster_count; i++)
        {
        double rg2 = 0.0;

        for (unsigned int j=0; j<m_clusters[i].size(); j++)
            {
            const unsigned int polymer_idx = m_clusters[i][j];

            for (unsigned int k=0; k<m_polymer_length; k++)
                {
                const unsigned int monomer_idx = h_polymer_index_map.data[m_polymer_indexer(k, polymer_idx)];
                Scalar3 dvec = make_scalar3(h_pos.data[monomer_idx].x - h_cluster_com.data[i].x,
                                            h_pos.data[monomer_idx].y - h_cluster_com.data[i].y,
                                            h_pos.data[monomer_idx].z - h_cluster_com.data[i].z);
                dvec = box.minImage(dvec);
                rg2 += dot(dvec, dvec);
                }
            }
        rg2 /= m_clusters[i].size();
        rg2 /= m_polymer_length;
        h_cluster_rg.data[i] = sqrt(rg2);
        }

    // Debug Output
    /*for (unsigned int i=0; i<m_polymer_count; i++)
        {
        std::cout<<"Polymer "<<i<<" is in cluster "<<m_cluster_ids[i]<<std::endl;
        }

    for (unsigned int i=0; i<m_cluster_count; i++)
        {
        std::cout<<"Cluster "<<i<<" contains "<<m_clusters[i].size()<<" polymers."<<std::endl;
        std::cout<<"Cluster center of mass is : "<<h_cluster_com.data[i].x<<", "<<h_cluster_com.data[i].y<<", "<<h_cluster_com.data[i].z<<std::endl;
        }*/
    }

/*! Function for calculating the forces between clusters. 
    Since the number of clusters will be comparatively low, use a simple for loop.
*/
void SurfaceCharge::Forces()
    {
    assert(m_pdata);
    assert(m_cluster_com);
    assert(m_cluster_rg);
    assert(m_polymer_index_map);
    
    // Acquire handle to required data and simulation box
    ArrayHandle<Scalar4> h_forces(m_force, access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_cluster_com(m_cluster_com, access_location::host, access_mode::read);
    ArrayHandle<Scalar> h_cluster_rg(m_cluster_rg, access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_polymer_index_map(m_polymer_index_map, access_location::host, access_mode::read);
    const BoxDim& box = m_pdata->getBox();

    vector<Scalar3> cluster_forces(m_cluster_count, make_scalar3(0.0, 0.0, 0.0));

    const Scalar pot_kappa2 = m_pot_kappa*m_pot_kappa;
    // Loop over all cluster pairs
    for (unsigned int i=0; i<m_cluster_count; i++)
        {
        for (unsigned int j=i+1; j<m_cluster_count; j++)
            {
            Scalar3 dvec = make_scalar3(h_cluster_com.data[j].x - h_cluster_com.data[i].x,
                                        h_cluster_com.data[j].y - h_cluster_com.data[i].y,
                                        h_cluster_com.data[j].z - h_cluster_com.data[i].z);
            dvec = box.minImage(dvec);

            const double sigma = h_cluster_rg.data[i] + h_cluster_rg.data[j];
            dvec.x /= sigma;
            dvec.y /= sigma;
            dvec.z /= sigma;
            const Scalar dr2 = dot(dvec, dvec);

            if (dr2 < m_pot_rcut2)
                {
                const double charge_i = 12.56637061*h_cluster_rg.data[i]*h_cluster_rg.data[i]; // prefactor is 4PI
                const double charge_j = 12.56637061*h_cluster_rg.data[j]*h_cluster_rg.data[j]; // prefactor is 4PI
                const double force = m_pot_epsilon*m_pot_epsilon*charge_i*charge_j/dr2/sigma;

                cluster_forces[i].x -= force*dvec.x;
                cluster_forces[i].y -= force*dvec.y;
                cluster_forces[i].z -= force*dvec.z;
                
                cluster_forces[j].x += force*dvec.x;
                cluster_forces[j].y += force*dvec.y;
                cluster_forces[j].z += force*dvec.z;
                }
            }
        }

    // First reset all forces of all particles
    for (unsigned int i=0; i<m_group->getNumMembers(); i++)
        {
        h_forces.data[i].x = 0.0;
        h_forces.data[i].y = 0.0;
        h_forces.data[i].z = 0.0;
        }

    // Now distribute forces in each cluster onto all monomers in that cluster
    for (unsigned int i=0; i<m_cluster_count; i++)
        {
        const Scalar3 force_fraction = make_scalar3(cluster_forces[i].x/m_clusters[i].size()/m_polymer_length,
                                                    cluster_forces[i].y/m_clusters[i].size()/m_polymer_length,
                                                    cluster_forces[i].z/m_clusters[i].size()/m_polymer_length);

        for (unsigned int j=0; j<m_clusters[i].size(); j++)
            {
            const unsigned int polymer_idx = m_clusters[i][j];

            //std::cout<<"Applying force ("<<force_fraction.x<<", "<<force_fraction.y<<", "<<force_fraction.z<<"), cluster:"<<i<<std::endl;

            for (unsigned int k=0; k<m_polymer_length; k++)
                {
                const unsigned int monomer_idx = h_polymer_index_map.data[m_polymer_indexer(k, polymer_idx)];
                h_forces.data[monomer_idx].x = force_fraction.x;
                h_forces.data[monomer_idx].y = force_fraction.y;
                h_forces.data[monomer_idx].z = force_fraction.z;

                //std::cout<<"Applying force ("<<force_fraction.x<<", "<<force_fraction.y<<", "<<force_fraction.z<<") on monomer: "<<monomer_idx<<" in polymer: "<<polymer_idx<<" in cluster: "<<i<<std::endl;
                }
            }
        }
    }

/*! Recursive function for assignig cluster ids
*/
void SurfaceCharge::MarkParticleAndNeighbors(unsigned int idx, unsigned int* cluster_id)
    {
    if (m_cluster_ids[idx] == -1)
        {
        m_cluster_ids[idx] = *cluster_id;

        for (unsigned int i=0; i<m_neighbors[idx].size(); i++)
            {
            MarkParticleAndNeighbors(m_neighbors[idx][i], cluster_id);
            }
        }
    }

/*! Perform the needed calculations to compute the forces
    \param timestep Current time step of the simulation
*/
void SurfaceCharge::computeForces(unsigned int timestep)
    {
    if (m_prof) m_prof->push("SurfaceCharge");

    // First check whether all parameters have been set
    if ((!m_cluster_params_set) || (!m_force_params_set))
        {
        std::cerr<<"SurfaceCharge parameters not set. Set using set_cluster_params and set_force_params."<<std::endl;
        return;
        }
    
    // If particles have been resorted, update mapping
    if (m_particles_sorted)
        {
        RemapPolymers();
        }

    // Calculate polymer center of masses
    CalcCenterOfMasses();

    // Perform cluster analysis
    ClusterAnalysis();

    // Calculate forces
    Forces();

    if (m_prof) m_prof->pop();
    }

/*! Function for setting cluster parameters
*/
void SurfaceCharge::setClusterParams(Scalar rcut)
    {
    m_cluster_rcut2 = rcut*rcut;
    m_cluster_params_set = true;
    std::cout<<"Cluster cutoff set to rcut="<<rcut<<std::endl;
    }

/*! Function for setting force parameters
*/
void SurfaceCharge::setForceParams(Scalar pot_epsilon, Scalar pot_kappa, Scalar pot_rcut)
    {
    m_pot_epsilon = pot_epsilon;
    m_pot_kappa = pot_kappa;
    m_pot_rcut2 = pot_rcut*pot_rcut;
    m_force_params_set = true;
    std::cout<<"Potential parameters set to: epsilon="<<m_pot_epsilon<<", kappa="<<m_pot_kappa<<", rcut="<<sqrt(m_pot_rcut2)<<std::endl;
    }

void export_SurfaceCharge()
    {
    class_<SurfaceCharge, boost::shared_ptr<SurfaceCharge>, bases<ForceCompute>, boost::noncopyable>
    ("SurfaceCharge", init< boost::shared_ptr<SystemDefinition>, boost::shared_ptr<ParticleGroup>, unsigned int >())
    .def("setClusterParams", &SurfaceCharge::setClusterParams)
    .def("setForceParams", &SurfaceCharge::setForceParams)
    ;
    }

