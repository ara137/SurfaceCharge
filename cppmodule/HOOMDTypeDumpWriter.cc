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

// Maintainer: joaander

/*! \file HOOMDTypeDumpWriter.cc
    \brief Defines the HOOMDTypeDumpWriter class
*/

#ifdef WIN32
#pragma warning( push )
#pragma warning( disable : 4244 )
#endif

#include <boost/python.hpp>
using namespace boost::python;

// begin mod
#include <vector>
// end mod

#include <sstream>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <boost/shared_ptr.hpp>

#include "HOOMDTypeDumpWriter.h"
#include "BondedGroupData.h"
#include "WallData.h"

#ifdef ENABLE_MPI
#include "Communicator.h"
#endif

using namespace std;
using namespace boost;

/*! \param sysdef SystemDefinition containing the ParticleData to dump
    \param base_fname The base name of the file xml file to output the information

    \note .timestep.xml will be apended to the end of \a base_fname when analyze() is called.
*/
HOOMDTypeDumpWriter::HOOMDTypeDumpWriter(boost::shared_ptr<SystemDefinition> sysdef, std::string base_fname)
        : Analyzer(sysdef), m_base_fname(base_fname), m_output_position(true),
        m_output_image(false), m_output_velocity(false), m_output_mass(false), m_output_diameter(false),
        m_output_type(false), m_output_bond(false), m_output_angle(false), m_output_wall(false),
        m_output_dihedral(false), m_output_improper(false), m_output_accel(false), m_output_body(false),
        m_output_charge(false), m_output_orientation(false), m_output_moment_inertia(false), m_vizsigma_set(false),
        writeMode(0), preserveAllTypes(true)
    {
    m_exec_conf->msg->notice(5) << "Constructing HOOMDTypeDumpWriter: " << base_fname << endl;
    }

HOOMDTypeDumpWriter::~HOOMDTypeDumpWriter()
    {
    m_exec_conf->msg->notice(5) << "Destroying HOOMDTypeDumpWriter" << endl;
    }

// begin mod
void HOOMDTypeDumpWriter::includeType(const std::string & type)
    {
    if(writeMode > -1)
        {
        assert(m_pdata);
        included.insert(m_pdata->getTypeByName(type));
        writeMode = 1;
        }
    else
        {
        m_exec_conf->msg->error() << "evaporator.xml: Cannot include AND exclude" << endl;
        throw runtime_error("Error choosing types for HOOMD dump file");
        }
        
    }

// void HOOMDTypeDumpWriter::includeTypes(const std::vector<string>& types)
//     {
//     if(writeMode > -1)
//         {
//         for(unsigned int i=0; i < types.size(); ++i)
//             {
//             assert(m_pdata);
//             included.insert(m_pdata->getTypeByName(types[i]));
//             }
//         writeMode = 1;
//         }
//     else
//         {
//         m_exec_conf->msg->error() << "evaporator.xml: Cannot include AND exclude" << endl;
//         throw runtime_error("Error choosing types for HOOMD dump file");
//         }
//     }

void HOOMDTypeDumpWriter::excludeType(const std::string & type)
    {
    if(writeMode < 1)
        {
        assert(m_pdata);
        excluded.insert(m_pdata->getTypeByName(type));
        writeMode = -1;
        }
    else
        {
        m_exec_conf->msg->error() << "evaporator.xml: Cannot include AND exclude" << endl;
        throw runtime_error("Error choosing types for HOOMD dump file");
        }
    }

// void HOOMDTypeDumpWriter::excludeTypes(const std::vector<string>& types)
//     {
//     if(writeMode < 1)
//         {
//         for(unsigned int i=0; i < types.size(); ++i)
//             {
//             assert(m_pdata);
//             excluded.insert(m_pdata->getTypeByName(types[i]));
//             }
//         writeMode = -1;
//         }
//     else
//         {
//         m_exec_conf->msg->error() << "evaporator.xml: Cannot include AND exclude" << endl;
//         throw runtime_error("Error choosing types for HOOMD dump file");
//         }
//     }

void HOOMDTypeDumpWriter::preserveTypes(bool enable)
    {
    preserveAllTypes = enable;
    }
//end mod

/*! \param enable Set to true to enable the writing of particle positions to the files in analyze()
*/
void HOOMDTypeDumpWriter::setOutputPosition(bool enable)
    {
    m_output_position = enable;
    }

/*! \param enable Set to true to enable the writing of particle images to the files in analyze()
*/
void HOOMDTypeDumpWriter::setOutputImage(bool enable)
    {
    m_output_image = enable;
    }

/*!\param enable Set to true to output particle velocities to the XML file on the next call to analyze()
*/
void HOOMDTypeDumpWriter::setOutputVelocity(bool enable)
    {
    m_output_velocity = enable;
    }

/*!\param enable Set to true to output particle masses to the XML file on the next call to analyze()
*/
void HOOMDTypeDumpWriter::setOutputMass(bool enable)
    {
    m_output_mass = enable;
    }

/*!\param enable Set to true to output particle diameters to the XML file on the next call to analyze()
*/
void HOOMDTypeDumpWriter::setOutputDiameter(bool enable)
    {
    m_output_diameter = enable;
    }

/*! \param enable Set to true to output particle types to the XML file on the next call to analyze()
*/
void HOOMDTypeDumpWriter::setOutputType(bool enable)
    {
    m_output_type = enable;
    }
/*! \param enable Set to true to output bonds to the XML file on the next call to analyze()
*/
void HOOMDTypeDumpWriter::setOutputBond(bool enable)
    {
    m_output_bond = enable;
    }
/*! \param enable Set to true to output angles to the XML file on the next call to analyze()
*/
void HOOMDTypeDumpWriter::setOutputAngle(bool enable)
    {
    m_output_angle = enable;
    }
/*! \param enable Set to true to output walls to the XML file on the next call to analyze()
*/
void HOOMDTypeDumpWriter::setOutputWall(bool enable)
    {
    m_output_wall = enable;
    }
/*! \param enable Set to true to output dihedrals to the XML file on the next call to analyze()
*/
void HOOMDTypeDumpWriter::setOutputDihedral(bool enable)
    {
    m_output_dihedral = enable;
    }
/*! \param enable Set to true to output impropers to the XML file on the next call to analyze()
*/
void HOOMDTypeDumpWriter::setOutputImproper(bool enable)
    {
    m_output_improper = enable;
    }
/*! \param enable Set to true to output acceleration to the XML file on the next call to analyze()
*/
void HOOMDTypeDumpWriter::setOutputAccel(bool enable)
    {
    m_output_accel = enable;
    }
/*! \param enable Set to true to output body to the XML file on the next call to analyze()
*/
void HOOMDTypeDumpWriter::setOutputBody(bool enable)
    {
    m_output_body = enable;
    }

/*! \param enable Set to true to output body to the XML file on the next call to analyze()
*/
void HOOMDTypeDumpWriter::setOutputCharge(bool enable)
    {
    m_output_charge = enable;
    }

/*! \param enable Set to true to output orientation to the XML file on the next call to analyze()
*/
void HOOMDTypeDumpWriter::setOutputOrientation(bool enable)
    {
    m_output_orientation = enable;
    }

/*! \param enable Set to true to output moment_inertia to the XML file on the next call to analyze()
*/
void HOOMDTypeDumpWriter::setOutputMomentInertia(bool enable)
    {
    m_output_moment_inertia = enable;
    }

/*! \param fname File name to write
    \param timestep Current time step of the simulation
*/
void HOOMDTypeDumpWriter::writeFile(std::string fname, unsigned int timestep)
    {
    // acquire the particle data
    SnapshotParticleData snapshot(m_pdata->getNGlobal());

    m_pdata->takeSnapshot(snapshot);

    BondData::Snapshot bdata_snapshot(m_sysdef->getBondData()->getNGlobal());
    if (m_output_bond) m_sysdef->getBondData()->takeSnapshot(bdata_snapshot);

    AngleData::Snapshot adata_snapshot(m_sysdef->getAngleData()->getNGlobal());
    if (m_output_angle) m_sysdef->getAngleData()->takeSnapshot(adata_snapshot);

    DihedralData::Snapshot ddata_snapshot(m_sysdef->getDihedralData()->getNGlobal());
    if (m_output_dihedral) m_sysdef->getDihedralData()->takeSnapshot(ddata_snapshot);

    ImproperData::Snapshot idata_snapshot(m_sysdef->getImproperData()->getNGlobal());
    if (m_output_improper) m_sysdef->getImproperData()->takeSnapshot(idata_snapshot);

#ifdef ENABLE_MPI
    // only the root processor writes the output file
    if (m_pdata->getDomainDecomposition() && ! m_exec_conf->isRoot())
        return;
#endif

    // begin mod
    std::vector<int> atomIdMap;
    std::vector<bool> typePreserved;
    unsigned int atomCounter = 0;
    atomIdMap.assign(m_pdata->getNGlobal(),-1);
    typePreserved.assign(m_pdata->getNTypes(),false);
    for(unsigned int j=0; j < m_pdata->getNGlobal(); ++j)
        {
        if(writeMode > 0) // included
            {
            if(included.find(snapshot.type[j]) != included.end() ||
              (preserveAllTypes && !typePreserved[snapshot.type[j]]) )
                {
                atomIdMap[j] = atomCounter;
                ++atomCounter;
                typePreserved[snapshot.type[j]] = true;
                }
            }
        else if(writeMode < 0) // excluded
            {
            if( excluded.find(snapshot.type[j]) == excluded.end() ||
                (preserveAllTypes && !typePreserved[snapshot.type[j]]) )
                {
                atomIdMap[j] = atomCounter;
                ++atomCounter;
                typePreserved[snapshot.type[j]] = true;
                }
            }
        else
            {
            atomIdMap[j] = atomCounter;
            ++atomCounter;
            typePreserved[snapshot.type[j]] = true;
            }
        }
    // end mod

    // open the file for writing
    ofstream f(fname.c_str());

    if (!f.good())
        {
        m_exec_conf->msg->error() << "dump.xml: Unable to open dump file for writing: " << fname << endl;
        throw runtime_error("Error writting hoomd_xml dump file");
        }

    BoxDim box = m_pdata->getGlobalBox();
    Scalar3 L = box.getL();
    Scalar xy = box.getTiltFactorXY();
    Scalar xz = box.getTiltFactorXZ();
    Scalar yz = box.getTiltFactorYZ();

    f.precision(13);
    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << "\n";
    f << "<hoomd_xml version=\"1.5\">" << "\n";
    f << "<configuration time_step=\"" << timestep << "\" "
      << "dimensions=\"" << m_sysdef->getNDimensions() << "\" "
    // begin mod
      << "natoms=\"" << atomCounter << "\" ";
    // end mod
    if (m_vizsigma_set)
        f << "vizsigma=\"" << m_vizsigma << "\" ";
    f << ">" << "\n";
    f << "<box " << "lx=\"" << L.x << "\" ly=\""<< L.y << "\" lz=\""<< L.z
      << "\" xy=\"" << xy << "\" xz=\"" << xz << "\" yz=\"" << yz << "\"/>" << "\n";

    f.precision(12);

    // If the position flag is true output the position of all particles to the file
    if (m_output_position)
        {
        // begin mod
        f << "<position num=\"" << atomCounter << "\">" << "\n";
        for (unsigned int j = 0; j < m_pdata->getNGlobal(); j++)
            {
            if(atomIdMap[j] > -1)
                {
                Scalar3 pos = snapshot.pos[j];

                f << pos.x << " " << pos.y << " "<< pos.z << "\n";

                if (!f.good())
                    {
                    m_exec_conf->msg->error() << "dump.xml: I/O error while writing HOOMD dump file" << endl;
                    throw runtime_error("Error writting HOOMD dump file");
                    }
                }
            }
        f <<"</position>" << "\n";
        // end mod
        }

    // If the image flag is true, output the image of each particle to the file
    if (m_output_image)
        {
        f << "<image num=\"" << atomCounter << "\">" << "\n";
        for (unsigned int j = 0; j < m_pdata->getNGlobal(); j++)
            {
            if(atomIdMap[j] > -1)
                {
                int3 image = snapshot.image[j];

                f << image.x << " " << image.y << " "<< image.z << "\n";

                if (!f.good())
                    {
                    m_exec_conf->msg->error() << "dump.xml: I/O error while writing HOOMD dump file" << endl;
                    throw runtime_error("Error writting HOOMD dump file");
                    }
                }
            }
        f <<"</image>" << "\n";
        }

    // If the velocity flag is true output the velocity of all particles to the file
    if (m_output_velocity)
        {
        f <<"<velocity num=\"" << atomCounter << "\">" << "\n";

        for (unsigned int j = 0; j < m_pdata->getNGlobal(); j++)
            {
            if(atomIdMap[j] > -1)
                {
                Scalar3 vel = snapshot.vel[j];
                f << vel.x << " " << vel.y << " " << vel.z << "\n";
                if (!f.good())
                    {
                    m_exec_conf->msg->error() << "dump.xml: I/O error while writing HOOMD dump file" << endl;
                    throw runtime_error("Error writting HOOMD dump file");
                    }
                }
            }

        f <<"</velocity>" << "\n";
        }

    // If the velocity flag is true output the velocity of all particles to the file
    if (m_output_accel)
        {
        f <<"<acceleration num=\"" << atomCounter << "\">" << "\n";

        for (unsigned int j = 0; j < m_pdata->getNGlobal(); j++)
            {
            if(atomIdMap[j] > -1)
                {
                Scalar3 accel = snapshot.accel[j];

                f << accel.x << " " << accel.y << " " << accel.z << "\n";
                if (!f.good())
                    {
                    m_exec_conf->msg->error() << "dump.xml: I/O error while writing HOOMD dump file" << endl;
                    throw runtime_error("Error writting HOOMD dump file");
                    }
                }
            }

        f <<"</acceleration>" << "\n";
        }

    // If the mass flag is true output the mass of all particles to the file
    if (m_output_mass)
        {
        f <<"<mass num=\"" << atomCounter << "\">" << "\n";

        for (unsigned int j = 0; j < m_pdata->getNGlobal(); j++)
            {
            if(atomIdMap[j] > -1)
                {
                Scalar mass = snapshot.mass[j];

                f << mass << "\n";
                if (!f.good())
                    {
                    m_exec_conf->msg->error() << "dump.xml: I/O error while writing HOOMD dump file" << endl;
                    throw runtime_error("Error writting HOOMD dump file");
                    }
                }
            }

        f <<"</mass>" << "\n";
        }

    // If the diameter flag is true output the mass of all particles to the file
    if (m_output_diameter)
        {
        f <<"<diameter num=\"" << atomCounter << "\">" << "\n";

        for (unsigned int j = 0; j < m_pdata->getNGlobal(); j++)
            {
            if(atomIdMap[j] > -1)
                {
                Scalar diameter = snapshot.diameter[j];
                f << diameter << "\n";
                if (!f.good())
                    {
                    m_exec_conf->msg->error() << "dump.xml: I/O error while writing HOOMD dump file" << endl;
                    throw runtime_error("Error writting HOOMD dump file");
                    }
                }
            }

        f <<"</diameter>" << "\n";
        }

    // If the Type flag is true output the types of all particles to an xml file
    if  (m_output_type)
        {
        f <<"<type num=\"" << atomCounter << "\">" << "\n";
        for (unsigned int j = 0; j < m_pdata->getNGlobal(); j++)
            {
            if(atomIdMap[j] > -1)
                {
                unsigned int type = snapshot.type[j];
                f << m_pdata->getNameByType(type) << "\n";
                }
            }
        f <<"</type>" << "\n";
        }

    // If the body flag is true output the bodies of all particles to an xml file
    // RIGID BODIES ARE NOT SUPPORTED YET
//     if  (m_output_body)
//         {
//         f <<"<body num=\"" << m_pdata->getNGlobal() << "\">" << "\n";
//         for (unsigned int j = 0; j < m_pdata->getNGlobal(); j++)
//             {
//             unsigned int body;
//             int out;
//             body = snapshot.body[j];
//             if (body == NO_BODY)
//                 out = -1;
//             else
//                 out = (int)body;
// 
//             f << out << "\n";
//             }
//         f <<"</body>" << "\n";
//         }

    // if the bond flag is true, output the bonds to the xml file
    if (m_output_bond)
        {
        f << "<bond num=\"" << bdata_snapshot.groups.size() << "\">" << "\n";
        boost::shared_ptr<BondData> bond_data = m_sysdef->getBondData();

        // loop over all bonds and write them out
        for (unsigned int i = 0; i < bdata_snapshot.groups.size(); i++)
            {
            BondData::members_t bond = bdata_snapshot.groups[i];
            int newTag0 = atomIdMap[bond.tag[0]];
            int newTag1 = atomIdMap[bond.tag[1]];
            if(newTag0 > -1 && newTag1 > -1)
                {
                unsigned int bond_type = bdata_snapshot.type_id[i];
                f << bond_data->getNameByType(bond_type) << " " << newTag0 << " " << newTag1 << "\n";
                }
            }

        f << "</bond>" << "\n";
        }

    // if the angle flag is true, output the angles to the xml file
    if (m_output_angle)
        {
        f << "<angle num=\"" << adata_snapshot.groups.size() << "\">" << "\n";
        boost::shared_ptr<AngleData> angle_data = m_sysdef->getAngleData();

        // loop over all angles and write them out
        for (unsigned int i = 0; i < adata_snapshot.groups.size(); i++)
            {
            AngleData::members_t angle = adata_snapshot.groups[i];
            int newTag0 = atomIdMap[angle.tag[0]];
            int newTag1 = atomIdMap[angle.tag[1]];
            int newTag2 = atomIdMap[angle.tag[2]];
            if(newTag0 > -1 && newTag1 > -1 && newTag2 > -1)
                {
                unsigned int angle_type = adata_snapshot.type_id[i];
                f << angle_data->getNameByType(angle_type) << " " << newTag0  << " " << newTag1 << " " << newTag2 << "\n";
                }
            }

        f << "</angle>" << "\n";
        }

    // if dihedral is true, write out dihedrals to the xml file
    if (m_output_dihedral)
        {
        f << "<dihedral num=\"" << ddata_snapshot.groups.size() << "\">" << "\n";
        boost::shared_ptr<DihedralData> dihedral_data = m_sysdef->getDihedralData();

        // loop over all angles and write them out
        for (unsigned int i = 0; i < ddata_snapshot.groups.size(); i++)
            {
            DihedralData::members_t dihedral = ddata_snapshot.groups[i];
            int newTag0 = atomIdMap[dihedral.tag[0]];
            int newTag1 = atomIdMap[dihedral.tag[1]];
            int newTag2 = atomIdMap[dihedral.tag[2]];
            int newTag3 = atomIdMap[dihedral.tag[3]];
            if(newTag0 > -1 && newTag1 > -1 && newTag2 > -1 && newTag3 > -1)
                {
                unsigned int dihedral_type = ddata_snapshot.type_id[i];
                f << dihedral_data->getNameByType(dihedral_type) << " " << newTag0  << " " << newTag1 << " "
                << newTag2 << " " << newTag3 << "\n";
                }
            }

        f << "</dihedral>" << "\n";
        }

    // if improper is true, write out impropers to the xml file
    if (m_output_improper)
        {
        f << "<improper num=\"" << idata_snapshot.groups.size() << "\">" << "\n";
        boost::shared_ptr<ImproperData> improper_data = m_sysdef->getImproperData();

        // loop over all angles and write them out
        for (unsigned int i = 0; i < idata_snapshot.groups.size(); i++)
            {
            ImproperData::members_t improper = idata_snapshot.groups[i];
            int newTag0 = atomIdMap[improper.tag[0]];
            int newTag1 = atomIdMap[improper.tag[1]];
            int newTag2 = atomIdMap[improper.tag[2]];
            int newTag3 = atomIdMap[improper.tag[3]];
            if(newTag0 > -1 && newTag1 > -1 && newTag2 > -1 && newTag3 > -1)
                {
                unsigned int improper_type = idata_snapshot.type_id[i];
                f << improper_data->getNameByType(improper_type) << " " << newTag0  << " " << newTag1 << " "
                << newTag2 << " " << newTag3 << "\n";
                }
            }

        f << "</improper>" << "\n";
        }

    // if the wall flag is true, output the walls to the xml file
    if (m_output_wall)
        {
        f << "<wall>" << "\n";
        boost::shared_ptr<WallData> wall_data = m_sysdef->getWallData();

        // loop over all walls and write them out
        for (unsigned int i = 0; i < wall_data->getNumWalls(); i++)
            {
            Wall wall = wall_data->getWall(i);
            f << "<coord ox=\"" << wall.origin_x << "\" oy=\"" << wall.origin_y << "\" oz=\"" << wall.origin_z <<
            "\" nx=\"" << wall.normal_x << "\" ny=\"" << wall.normal_y << "\" nz=\"" << wall.normal_z << "\" />" << "\n";
            }
        f << "</wall>" << "\n";
        }

    // If the charge flag is true output the mass of all particles to the file
    if (m_output_charge)
        {
        f <<"<charge num=\"" << atomCounter << "\">" << "\n";

        for (unsigned int j = 0; j < m_pdata->getNGlobal(); j++)
            {
            if(atomIdMap[j] > -1)
                {
                Scalar charge = snapshot.charge[j];
                f << charge << "\n";
                if (!f.good())
                    {
                    m_exec_conf->msg->error() << "dump.xml: I/O error while writing HOOMD dump file" << endl;
                    throw runtime_error("Error writting HOOMD dump file");
                    }
                }
            }

        f <<"</charge>" << "\n";
        }

    // if the orientation flag is set, write out the orientation quaternion to the XML file
    if (m_output_orientation)
        {
        f << "<orientation num=\"" << atomCounter << "\">" << "\n";

        for (unsigned int j = 0; j < m_pdata->getNGlobal(); j++)
            {
            if(atomIdMap[j] > -1)
                {
                // use the rtag data to output the particles in the order they were read in
                Scalar4 orientation = snapshot.orientation[j];
                f << orientation.x << " " << orientation.y << " " << orientation.z << " " << orientation.w << "\n";
                if (!f.good())
                    {
                    m_exec_conf->msg->error() << "dump.xml: I/O error while writing HOOMD dump file" << endl;
                    throw runtime_error("Error writting HOOMD dump file");
                    }
                }
            }
        f << "</orientation>" << "\n";
        }

    // if the moment_inertia flag is set, write out the orientation quaternion to the XML file
    if (m_output_moment_inertia)
        {
        f << "<moment_inertia num=\"" << atomCounter << "\">" << "\n";

        for (unsigned int i = 0; i < m_pdata->getNGlobal(); i++)
            {
            if(atomIdMap[i] > -1)
                {
                // inertia tensors are stored by tag
                InertiaTensor I = snapshot.inertia_tensor[i];
                for (unsigned int c = 0; c < 5; c++)
                    f << I.components[c] << " ";
                f << I.components[5] << "\n";

                if (!f.good())
                    {
                    m_exec_conf->msg->error() << "dump.xml: I/O error while writing HOOMD dump file" << endl;
                    throw runtime_error("Error writting HOOMD dump file");
                    }
                }
            }
        f << "</moment_inertia>" << "\n";
        }

    f << "</configuration>" << "\n";
    f << "</hoomd_xml>" << "\n";

    if (!f.good())
        {
        m_exec_conf->msg->error() << "dump.xml: I/O error while writing HOOMD dump file" << endl;
        throw runtime_error("Error writting HOOMD dump file");
        }

    f.close();

    }

/*! \param timestep Current time step of the simulation
    Writes a snapshot of the current state of the ParticleData to a hoomd_xml file.
*/
void HOOMDTypeDumpWriter::analyze(unsigned int timestep)
    {
    if (m_prof)
        m_prof->push("Dump XML");

    ostringstream full_fname;
    string filetype = ".xml";

    // Generate a filename with the timestep padded to ten zeros
    full_fname << m_base_fname << "." << setfill('0') << setw(10) << timestep << filetype;
    writeFile(full_fname.str(), timestep);

    if (m_prof)
        m_prof->pop();
    }

void export_HOOMDTypeDumpWriter()
    {   
    class_<HOOMDTypeDumpWriter, boost::shared_ptr<HOOMDTypeDumpWriter>, bases<Analyzer>, boost::noncopyable>
    ("HOOMDTypeDumpWriter", init< boost::shared_ptr<SystemDefinition>, std::string >())
    .def("includeType", &HOOMDTypeDumpWriter::includeType)
    .def("excludeType", &HOOMDTypeDumpWriter::excludeType)
    .def("preserveTypes", &HOOMDTypeDumpWriter::preserveTypes)
    .def("setOutputPosition", &HOOMDTypeDumpWriter::setOutputPosition)
    .def("setOutputImage", &HOOMDTypeDumpWriter::setOutputImage)
    .def("setOutputVelocity", &HOOMDTypeDumpWriter::setOutputVelocity)
    .def("setOutputMass", &HOOMDTypeDumpWriter::setOutputMass)
    .def("setOutputDiameter", &HOOMDTypeDumpWriter::setOutputDiameter)
    .def("setOutputType", &HOOMDTypeDumpWriter::setOutputType)
    .def("setOutputBody", &HOOMDTypeDumpWriter::setOutputBody)
    .def("setOutputBond", &HOOMDTypeDumpWriter::setOutputBond)
    .def("setOutputAngle", &HOOMDTypeDumpWriter::setOutputAngle)
    .def("setOutputDihedral", &HOOMDTypeDumpWriter::setOutputDihedral)
    .def("setOutputImproper", &HOOMDTypeDumpWriter::setOutputImproper)
    .def("setOutputWall", &HOOMDTypeDumpWriter::setOutputWall)
    .def("setOutputAccel", &HOOMDTypeDumpWriter::setOutputAccel)
    .def("setOutputCharge", &HOOMDTypeDumpWriter::setOutputCharge)
    .def("setOutputOrientation", &HOOMDTypeDumpWriter::setOutputOrientation)
    .def("setVizSigma", &HOOMDTypeDumpWriter::setVizSigma)
    .def("writeFile", &HOOMDTypeDumpWriter::writeFile)
    ;
    }

#ifdef WIN32
#pragma warning( pop )
#endif
