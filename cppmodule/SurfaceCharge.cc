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
SurfaceCharge::SurfaceCharge(boost::shared_ptr<SystemDefinition> sysdef, boost::shared_ptr<ParticleGroup> group)
        : ForceCompute(sysdef)
    {
    // assign group and resize polymer map
    m_group = group;
    m_polymer_index_map.resize(group->getNumMembers());
    // This should be changed later to a member variable which can be changed externally
    m_polymer_length = 32;
    
    // Check for stray monomers
    if (m_group->getNumMembers()%m_polymer_length != 0)
        {
        std::cout<<"SurfaceCharge: Number of particles is not an integer multiple of monomers/chain!"<<std::endl;
        }

    // Set polymer count
    m_polymer_count = m_group->getNumMembers()/m_polymer_length;
    // Resize center of masses array
    m_polymer_com.resize(m_polymer_count);
// Initialize indexer, where each row consists of one polymer
    m_polymer_indexer = Index2D(m_polymer_length, m_polymer_count);

    // Get simulation box for periodic boundary conditions
    const BoxDim& box = m_pdata->getGlobalBox();
    m_box_size = box.getL();
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

    std::cout<<"SurfaceCharge: Remapping Polymers"<<std::endl;
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
    ArrayHandle<unsigned int> h_polymer_index_map(m_polymer_index_map, access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_polymer_com(m_polymer_com, access_location::host, access_mode::readwrite);
    

    for (unsigned int i=0; i<m_polymer_count; i++)
        {
        Scalar4 com = make_scalar4(0.0, 0.0, 0.0, 0.0);
        // index of the previous bead
        unsigned int idx_prev; 

        for (unsigned int j=0; j<m_polymer_length; j++)
            {
            const unsigned int idx = m_polymer_indexer(i, j);
            com.x += h_pos.data[idx].x;
            com.y += h_pos.data[idx].y;
            com.z += h_pos.data[idx].z;

            // Take care of periodic boundary conditions
            if (j > 0)
                {
                if ((h_pos.data[idx].x - h_pos.data[idx_prev].x) > 0.5*m_box_size.x)
                    {
                    com.x -= m_box_size.x;
                    }
                else if ((h_pos.data[idx].x - h_pos.data[idx_prev].x) < -0.5*m_box_size.x)
                    {
                    com.x += m_box_size.x;
                    }
                }
                
            idx_prev = idx;
            }
        }
    }

/*! Perform the needed calculations to compute the forces
    \param timestep Current time step of the simulation
*/
void SurfaceCharge::computeForces(unsigned int timestep)
    {
    if (m_prof) m_prof->push("ExampleUpdater");
    
    // If particles have been resorted, update mapping
    if (m_particles_sorted)
        {
        RemapPolymers();
        }

    // Calculate polymer center of masses
    CalcCenterOfMasses();

    // access the particle data for writing on the CPU
    assert(m_pdata);
    assert(m_polymer_index_map);
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_polymer_index_map(m_polymer_index_map, access_location::host, access_mode::read);
    // check if particles have been sorted, and then update internal polymer map if necessary

    for (unsigned int i = 0; i < m_pdata->getN(); i++)
        {
        /*h_vel.data[i].x = Scalar(0.0);
        h_vel.data[i].y = Scalar(0.0);
        h_vel.data[i].z = Scalar(0.0);*/
        }

    if (m_prof) m_prof->pop();
    }

void export_SurfaceCharge()
    {
    class_<SurfaceCharge, boost::shared_ptr<SurfaceCharge>, bases<ForceCompute>, boost::noncopyable>
    ("SurfaceCharge", init< boost::shared_ptr<SystemDefinition>, boost::shared_ptr<ParticleGroup> >())
    ;
    }
