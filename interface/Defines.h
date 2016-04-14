#pragma once

// Debugging macros

// #define _TT_DEBUG_
#define TT_MTT_DEBUG (false)
#define TT_HLT_DEBUG (false)
#define TT_GEN_DEBUG (false)


#if TT_GEN_DEBUG
#define FILL_GEN_COLL( X ) \
    if (flags.isLastCopy()) { \
        genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
        gen_##X = gen_index; \
        gen_index++; \
        std::cout << "Assigning gen_" #X " = " << i << " (" << pdg_id << ")" << std::endl; \
    } \
    if (flags.isFirstCopy()) { \
        genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
        gen_##X##_beforeFSR = gen_index; \
        gen_index++; \
        std::cout << "Assigning gen_" #X "_beforeFSR = " << i << " (" << pdg_id << ")" << std::endl; \
    }
#else
#define FILL_GEN_COLL( X ) \
    if (flags.isLastCopy()) { \
        genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
        gen_##X = gen_index; \
        gen_index++; \
    } \
    if (flags.isFirstCopy()) { \
        genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
        gen_##X##_beforeFSR = gen_index; \
        gen_index++; \
    }
#endif

// Assign index to X if it's empty, or Y if not
#if TT_GEN_DEBUG
#define FILL_GEN_COLL2(X, Y, ERROR) \
    if (flags.isLastCopy()) { \
        if (gen_##X == -1){ \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##X = gen_index; \
            gen_index++; \
            std::cout << "Assigning gen_" #X " = " << i << " (" << pdg_id << ")" << std::endl; \
        } else if (gen_##Y == -1) { \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##Y = gen_index; \
            gen_index++; \
            std::cout << "Assigning gen_" #Y " = " << i << " (" << pdg_id << ")" << std::endl; \
        } else \
            std::cout << ERROR << std::endl; \
    } \
    if (flags.isFirstCopy()) { \
        if (gen_##X##_beforeFSR == -1) { \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##X##_beforeFSR = gen_index; \
            gen_index++; \
            std::cout << "Assigning gen_" #X "_beforeFSR = " << i << " (" << pdg_id << ")" << std::endl; \
        } else if (gen_##Y##_beforeFSR == -1) { \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##Y##_beforeFSR = gen_index; \
            gen_index++; \
            std::cout << "Assigning gen_" #Y "_beforeFSR = " << i << " (" << pdg_id << ")" << std::endl; \
        } else \
            std::cout << ERROR << std::endl; \
    }
#else
#define FILL_GEN_COLL2(X, Y, ERROR) \
    if (flags.isLastCopy()) { \
        if (gen_##X == -1){ \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##X = gen_index; \
            gen_index++; \
        } else if (gen_##Y == -1) { \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##Y = gen_index; \
            gen_index++; \
        } else \
            std::cout << ERROR << std::endl; \
    } \
    if (flags.isFirstCopy()) { \
        if (gen_##X##_beforeFSR == -1) { \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##X##_beforeFSR = gen_index; \
            gen_index++; \
        } else if (gen_##Y##_beforeFSR == -1) { \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##Y##_beforeFSR = gen_index; \
            gen_index++; \
        } else \
            std::cout << ERROR << std::endl; \
    }
#endif
