/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_MACROS_HPP
#define MCPACK_MACROS_HPP

#define MCPACK_MAJOR_VERSION 0
#define MCPACK_MINOR_VERSION 1

#ifndef MCPACK_NDEBUG
#   define MCPACK_ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::exit(EXIT_FAILURE); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif /*MCPACK_NDEBUG */

#endif //MCPACK_MACROS_HPP
