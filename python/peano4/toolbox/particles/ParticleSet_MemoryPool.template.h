//
// Peano4 data file
// Generated by Peano's Python API
// www.peano-framework.org
// This is generated. Be careful with adding your own stuff
//
#pragma once


#include <list>

#include "tarch/logging/Log.h"
#include "tarch/la/Vector.h"
#include "tarch/mpi/mpi.h"

#include "peano4/utils/Globals.h"
#include "peano4/datamanagement/VertexMarker.h"
#include "peano4/grid/TraversalObserver.h"

#include "toolbox/particles/ParticleSet.h"
#include "toolbox/particles/memorypool/VertexWiseContinuousMemoryPool.h"
#include "toolbox/particles/memorypool/GlobalContinuousMemoryPool.h"


{% for item in NAMESPACE -%}
  {% if not loop.last %}
  namespace {{ item }} {
  {% endif %}
{%- endfor %}

  namespace globaldata {
    class {{PARTICLE_TYPE}};
  }
{% for item in NAMESPACE -%}
  {% if not loop.last %}
  }
  {% endif %}
{%- endfor %}



{% for item in NAMESPACE -%}
  namespace {{ item }} {
{%- endfor %}


  class {{CLASSNAME}}: public toolbox::particles::ParticleSet< {% for item in NAMESPACE -%}{% if not loop.last %}{{ item }}::{% endif %}{%- endfor %}globaldata::{{PARTICLE_TYPE}} >
                       {
    public:
      using Base       = toolbox::particles::ParticleSet< {% for item in NAMESPACE -%}{% if not loop.last %}{{ item }}::{% endif %}{%- endfor %}globaldata::{{PARTICLE_TYPE}} >;
      using MemoryPool = {{MEMORY_POOL_TYPE}}< {% for item in NAMESPACE -%}{% if not loop.last %}{{ item }}::{% endif %}{%- endfor %}globaldata::{{PARTICLE_TYPE}} >;

      /**
       * This is the container as it looks like from outside.
       * Internally, we use instances of MemoryPool to administer
       * the data.
       */
      typedef MemoryPool::Container          Container;

      /**
       * Expose C++ standard interface
       */
      typedef typename Container::value_type      value_type;
      typedef typename Container::iterator        iterator;
      typedef typename Container::const_iterator  const_iterator;

      {{CLASSNAME}}() = default;


      enum ObjectConstruction {
        NoData
      };

      {{CLASSNAME}}( ObjectConstruction ):
        {{CLASSNAME}}() {}

      /**
       * Shallow copy constructor
       *
       * There's no need to redefine this one. I just wanted to highlight
       * that it is a shallow copy, i.e. to add a comment, and so I
       * "overwrote" it.
       *
       * Please consult @ref toolbox_particles_memorypool "the generic discussion of the particle memory management" for
       * further context.
       */
      {{CLASSNAME}}( const {{CLASSNAME}}& ) = default;

      /**
       * Merge two particle sets
       *
       * Please consult @ref toolbox_particles_memorypool "the generic discussion of the particle memory management" before
       * you study the implementation remarks below. Furhermore, consult the
       * documentation of peano4.toolbox.particles.api.UpdateParallelState which
       * provides some examples of the functionality added here.
       *
       * @see STDVectorOverContainerOfPointers::startReceive()
       * @see clone()
       */
      void merge(::peano4::grid::TraversalObserver::SendReceiveContext context, {{CLASSNAME}}& neighbour, const peano4::datamanagement::VertexMarker& marker, int spacetreeId);

      static bool send(const peano4::datamanagement::VertexMarker& marker);
      static bool receiveAndMerge(const peano4::datamanagement::VertexMarker& marker);
      static ::peano4::grid::LoadStoreComputeFlag loadStoreComputeFlag(const peano4::datamanagement::VertexMarker& marker);

      /**
       * Deep clone
       *
       * Please consult @ref toolbox_particles_memorypool "the generic discussion of the particle memory management" before
       * you study the implementation remarks below.
       *
       * @see STDVectorOverContainerOfPointers() for a description of the
       *   intra-node parallelisation. This description provides an overview
       *   of the overall lifecycle of the data.
       * @see merge() for the corresponding delete that eliminates the deep
       *   copies.
       */
      void clone( const {{CLASSNAME}}& otherSet );

      std::string toString() const;

      /**
       * Mark whole particle set as set of particles which can't be lifted
       *
       * Consult @ref toolbox_particles_memorypool "the generic discussion of the particle memory management" as well as the documentation of toolbox::particles::ParticleSet for
       * some context and references to other routines.
       *
       * The operation runs over all hosted particles and invokes
       * particleCanNotBeLiftedLocally() per particle.
       */
      void particlesCanNotBeLiftedLocally();

      /**
       * Shallow clear
       *
       * Plain delegate to superclass clear(), i.e. this routine clears a set
       * of pointers, but it does not delete the actual particles. If you want
       * that behaviour, you have to call deleteParticles().
       *
       * As it is a shallow clear, we may assume that another stack holds a
       * copy to "our" particles. What can now happen is the following: A mesh
       * traversal consolidates all data structures, i.e. invokes gather() on
       * its vertices. At the end, these gathered vertices are also pushed onto
       * to parallel data exchange stacks. Here, we employ the default copy
       * constructor/assignment operator. That means, the boundary stack now
       * "thinks" it holds a continuous piece of memory. Actually, it does so,
       * but it does not own this memory. Therefore, the clear() clears indeed
       * the list of pointers, and it clears the container, but it does not
       * free the underlying memory.
       */
      void clear();

      /**
       * Inform set that particle cannot be lifted on its current tree
       *
       * Consult @ref toolbox_particles_memorypool "the generic discussion of the particle memory management" as well as the documentation of toolbox::particles::ParticleSet for
       * some context and references to other routines.
       */
      Container::iterator particleCanNotBeLiftedLocally( const Container::iterator&   particle );

      void particlesCanNotBeDroppedLocallySoRollOverForNextMeshSweep(DoFType*  particle);

      /**
       * This is a deep delete
       *
       * Consult @ref toolbox_particles_memorypool "the generic discussion of the particle memory management" as well as the documentation of toolbox::particles::ParticleSet for
       * some context and references to other routines.
       */
      void deleteParticles();

      /**
       * Expose the STD container partially
       *
       * The fact that you hold a const or non-const iterator does not mean that
       * you are ever allowed to delete a particle. That has to be done through
       * the corresponding deleteParticle() operation.
       *
       * Consult @ref toolbox_particles_memorypool "the generic discussion of the particle memory management" as well as the documentation of toolbox::particles::ParticleSet for
       * some context and references to other routines.
       */
      Container::iterator begin();
      Container::iterator end();
      Container::const_iterator begin() const;
      Container::const_iterator end() const;
      int size() const;
      bool empty() const;

      /**
       * Add one particle
       *
       * Consult @ref toolbox_particles_memorypool "the generic discussion of the particle memory management" as well as the documentation of toolbox::particles::ParticleSet for
       * some context and references to other routines.
       *
       * This routine does not copy the particle. It adds the references
       * (pointer) to the particle only. Notably, this routine does not
       * give the set exclusive ownership of the particle.
       *
       * Is internally mapped onto a push_back.
       *
       * @param eliminateDuplicates If particleSet contains a duplicate that's
       *   already assigned to the vertex, ignore this duplicate.
       */
      void addParticle(DoFType* particle);

      /**
       * Some simple validation
       *
       * I run some simple validations such as: does every particle have a
       * valid parallel state. This might result in the detection of memory
       * garbage. The routine is typically used within assertions.
       */
      bool isValid() const;

      /**
       * Moves particles from one particle set to the other
       *
       * This is used by the multiscale algorithm
       */
      void grabParticles({{CLASSNAME}}& sourceParticleSet);

      /**
       * Move particle from one container to the other one
       *
       * ## Difference to particles on heap
       *
       * I need a copy of the particle iterator, as I might have to scatter
       * sourceParticleSet and then redirect particle to an entry within a
       * scattered set.
       *
       * @return Iterator within sourceParticleSet which is returned after
       *   erasing it from this set.
       */
      Container::iterator grabParticle(
        Container::iterator  particle,
        {{CLASSNAME}}&       sourceParticleSet
      );

      /**
       * Delete a particle
       *
       * Should be used, for example, if a particle has left the local domain.
       * Consult @ref toolbox_particles_memorypool "the generic discussion of the particle memory management" as well as the documentation of toolbox::particles::ParticleSet for
       * some context and references to other routines.
       * 
       *
       * ## Difference to particles on heap
       *
       * I need a copy of the particle iterator, as I might have to scatter
       * sourceParticleSet and then redirect particle to an entry within a
       * scattered set.
       *
       *
       * @return an iterator to the particle. This iterator is robust, i.e. even 
       *   if we have altered the underlying data, this one is still valid. It is
       *   however not efficient, i.e. it might point to begin() and any outer
       *   loop hence might restart to run over all particles again. 
       */
      Container::iterator deleteParticle( Container::iterator  particle );

      /**
       * Gather particles
       *
       * This routine consolidates all the associated particles into one
       * chunk. This is one of the two routines that we offer compared to a
       * realisation which scatters particles over the heap.
       */
      void gather();

      /**
       * Scatter data
       *
       * This routine is usually not called on the user side, as we try to
       * hide it, i.e. invoke scatter internally. There are very few cases
       * where you might want to scatter manually.
       */
      void scatter();

      /**
       * Are particles gathered atm
       *
       * Just a query. This is one of the two routines that we offer compared to a
       * realisation which scatters particles over the heap.
       */
      bool isGathered() const;

      /**
       * Extracted into routine of its own
       *
       * This routine it used by boundary merges, periodic boundary merges, and
       * also drops.
       *
       * ## Debug data
       *
       * Even if we replace a current particle with an incoming one, we keep
       * the current particle's (debug) location. This way, we ensure that all
       * debug data remain correct also for periodic BCs.
       *
       * @param neighbourParticle Will either be incorporated into current data
       *   structure or the routine will delete it if it decides not to use
       *   the object the pointer points to.
       */
      void mergeWithParticle(
        DoFType*                                     neighbourParticle,
        const peano4::datamanagement::VertexMarker&  marker,
        int                                          spacetreeId
      );
    private:
      static tarch::logging::Log _log;

      MemoryPool  _memoryPool;
};

{% for item in NAMESPACE -%}
}
{%- endfor %}


