##
##  Copyright © 2009 CNRS
##  Copyright © 2009-2018 Inria.  All rights reserved.
##  Copyright © 2009-2012 Université Bordeaux
##  Copyright © 2009-2011 Cisco Systems, Inc.  All rights reserved.
##  See COPYING in top-level directory.
##
## * \file
##  \brief The bitmap API, for use in hwloc itself.
##

## * \defgroup hwlocality_bitmap The bitmap API
##
##  The ::hwloc_bitmap_t type represents a set of integers (positive or null).
##  A bitmap may be of infinite size (all bits are set after some point).
##  A bitmap may even be full if all bits are set.
##
##  Bitmaps are used by hwloc for sets of OS processors
##  (which may actually be hardware threads) as by ::hwloc_cpuset_t
##  (a typedef for ::hwloc_bitmap_t), or sets of NUMA memory nodes
##  as ::hwloc_nodeset_t (also a typedef for ::hwloc_bitmap_t).
##  Those are used for cpuset and nodeset fields in the ::hwloc_obj structure,
##  see \ref hwlocality_object_sets.
##
##  <em>Both CPU and node sets are always indexed by OS physical number.</em>
##  However users should usually not build CPU and node sets manually
##  (e.g. with hwloc_bitmap_set()).
##  One should rather use existing object sets and combine them with
##  hwloc_bitmap_or(), etc.
##  For instance, binding the current thread on a pair of cores may be performed with:
##  \code
##  hwloc_obj_t core1 = ... , core2 = ... ;
##  hwloc_bitmap_t set = hwloc_bitmap_alloc();
##  hwloc_bitmap_or(set, core1->cpuset, core2->cpuset);
##  hwloc_set_cpubind(topology, set, HWLOC_CPUBIND_THREAD);
##  hwloc_bitmap_free(set);
##  \endcode
##
##  \note Most functions below return an int that may be negative in case of
##  error. The usual error case would be an internal failure to realloc/extend
##  the storage of the bitmap (\p errno would be set to \c ENOMEM).
##
##  \note Several examples of using the bitmap API are available under the
##  doc/examples/ directory in the source tree.
##  Regression tests such as tests/hwloc/hwloc_bitmap*.c also make intensive use
##  of this API.
##  @{
##

when hostOs == "macosx":
  const hwloclib {.strDefine.} = "libhwloc.dylib"
else:
  const hwloclib {.strDefine.} = "libhwloc.so"
{.pragma: hwloc, dynlib:hwloclib.}

type
  hwloc_bitmap_s {.bycopy.} = object

## * \brief
##  Set of bits represented as an opaque pointer to an internal bitmap.
##

type
  hwloc_bitmap_t* = ptr hwloc_bitmap_s

## * \brief a non-modifiable ::hwloc_bitmap_t

type
  hwloc_const_bitmap_t* = ptr hwloc_bitmap_s

##
##  Bitmap allocation, freeing and copying.
##
## * \brief Allocate a new empty bitmap.
##
##  \returns A valid bitmap or \c NULL.
##
##  The bitmap should be freed by a corresponding call to
##  hwloc_bitmap_free().
##

proc hwloc_bitmap_alloc*(): hwloc_bitmap_t {.importc: "hwloc_bitmap_alloc",
    hwloc.}
## * \brief Allocate a new full bitmap.

proc hwloc_bitmap_alloc_full*(): hwloc_bitmap_t {.
    importc: "hwloc_bitmap_alloc_full", hwloc.}
## * \brief Free bitmap \p bitmap.
##
##  If \p bitmap is \c NULL, no operation is performed.
##

proc hwloc_bitmap_free*(bitmap: hwloc_bitmap_t) {.importc: "hwloc_bitmap_free",
    hwloc.}
## * \brief Duplicate bitmap \p bitmap by allocating a new bitmap and copying \p bitmap contents.
##
##  If \p bitmap is \c NULL, \c NULL is returned.
##

proc hwloc_bitmap_dup*(bitmap: hwloc_const_bitmap_t): hwloc_bitmap_t {.
    importc: "hwloc_bitmap_dup", hwloc.}
## * \brief Copy the contents of bitmap \p src into the already allocated bitmap \p dst

proc hwloc_bitmap_copy*(dst: hwloc_bitmap_t; src: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_copy", hwloc.}
##
##  Bitmap/String Conversion
##
## * \brief Stringify a bitmap.
##
##  Up to \p buflen characters may be written in buffer \p buf.
##
##  If \p buflen is 0, \p buf may safely be \c NULL.
##
##  \return the number of character that were actually written if not truncating,
##  or that would have been written (not including the ending \\0).
##

proc hwloc_bitmap_snprintf*(buf: cstring; buflen: csize_t; bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_snprintf", hwloc.}
## * \brief Stringify a bitmap into a newly allocated string.
##
##  \return -1 on error.
##

proc hwloc_bitmap_asprintf*(strp: ptr cstring; bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_asprintf", hwloc.}
## * \brief Parse a bitmap string and stores it in bitmap \p bitmap.
##

proc hwloc_bitmap_sscanf*(bitmap: hwloc_bitmap_t; string: cstring): cint {.
    importc: "hwloc_bitmap_sscanf", hwloc.}
## * \brief Stringify a bitmap in the list format.
##
##  Lists are comma-separated indexes or ranges.
##  Ranges are dash separated indexes.
##  The last range may not have an ending indexes if the bitmap is infinitely set.
##
##  Up to \p buflen characters may be written in buffer \p buf.
##
##  If \p buflen is 0, \p buf may safely be \c NULL.
##
##  \return the number of character that were actually written if not truncating,
##  or that would have been written (not including the ending \\0).
##

proc hwloc_bitmap_list_snprintf*(buf: cstring; buflen: csize_t;
                                bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_list_snprintf", hwloc.}
## * \brief Stringify a bitmap into a newly allocated list string.
##
##  \return -1 on error.
##

proc hwloc_bitmap_list_asprintf*(strp: ptr cstring; bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_list_asprintf", hwloc.}
## * \brief Parse a list string and stores it in bitmap \p bitmap.
##

proc hwloc_bitmap_list_sscanf*(bitmap: hwloc_bitmap_t; string: cstring): cint {.
    importc: "hwloc_bitmap_list_sscanf", hwloc.}
## * \brief Stringify a bitmap in the taskset-specific format.
##
##  The taskset command manipulates bitmap strings that contain a single
##  (possible very long) hexadecimal number starting with 0x.
##
##  Up to \p buflen characters may be written in buffer \p buf.
##
##  If \p buflen is 0, \p buf may safely be \c NULL.
##
##  \return the number of character that were actually written if not truncating,
##  or that would have been written (not including the ending \\0).
##

proc hwloc_bitmap_taskset_snprintf*(buf: cstring; buflen: csize_t;
                                   bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_taskset_snprintf", hwloc.}
## * \brief Stringify a bitmap into a newly allocated taskset-specific string.
##
##  \return -1 on error.
##

proc hwloc_bitmap_taskset_asprintf*(strp: ptr cstring;
                                   bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_taskset_asprintf", hwloc.}
## * \brief Parse a taskset-specific bitmap string and stores it in bitmap \p bitmap.
##

proc hwloc_bitmap_taskset_sscanf*(bitmap: hwloc_bitmap_t; string: cstring): cint {.
    importc: "hwloc_bitmap_taskset_sscanf", hwloc.}
##
##  Building bitmaps.
##
## * \brief Empty the bitmap \p bitmap

proc hwloc_bitmap_zero*(bitmap: hwloc_bitmap_t) {.importc: "hwloc_bitmap_zero",
    hwloc.}
## * \brief Fill bitmap \p bitmap with all possible indexes (even if those objects don't exist or are otherwise unavailable)

proc hwloc_bitmap_fill*(bitmap: hwloc_bitmap_t) {.importc: "hwloc_bitmap_fill",
    hwloc.}
## * \brief Empty the bitmap \p bitmap and add bit \p id

proc hwloc_bitmap_only*(bitmap: hwloc_bitmap_t; id: cuint): cint {.
    importc: "hwloc_bitmap_only", hwloc.}
## * \brief Fill the bitmap \p and clear the index \p id

proc hwloc_bitmap_allbut*(bitmap: hwloc_bitmap_t; id: cuint): cint {.
    importc: "hwloc_bitmap_allbut", hwloc.}
## * \brief Setup bitmap \p bitmap from unsigned long \p mask

proc hwloc_bitmap_from_ulong*(bitmap: hwloc_bitmap_t; mask: culong): cint {.
    importc: "hwloc_bitmap_from_ulong", hwloc.}
## * \brief Setup bitmap \p bitmap from unsigned long \p mask used as \p i -th subset

proc hwloc_bitmap_from_ith_ulong*(bitmap: hwloc_bitmap_t; i: cuint; mask: culong): cint {.
    importc: "hwloc_bitmap_from_ith_ulong", hwloc.}
## * \brief Setup bitmap \p bitmap from unsigned longs \p masks used as first \p nr subsets

proc hwloc_bitmap_from_ulongs*(bitmap: hwloc_bitmap_t; nr: cuint; masks: ptr culong): cint {.
    importc: "hwloc_bitmap_from_ulongs", hwloc.}
##
##  Modifying bitmaps.
##
## * \brief Add index \p id in bitmap \p bitmap

proc hwloc_bitmap_set*(bitmap: hwloc_bitmap_t; id: cuint): cint {.
    importc: "hwloc_bitmap_set", hwloc.}
## * \brief Add indexes from \p begin to \p end in bitmap \p bitmap.
##
##  If \p end is \c -1, the range is infinite.
##

proc hwloc_bitmap_set_range*(bitmap: hwloc_bitmap_t; begin: cuint; `end`: cint): cint {.
    importc: "hwloc_bitmap_set_range", hwloc.}
## * \brief Replace \p i -th subset of bitmap \p bitmap with unsigned long \p mask

proc hwloc_bitmap_set_ith_ulong*(bitmap: hwloc_bitmap_t; i: cuint; mask: culong): cint {.
    importc: "hwloc_bitmap_set_ith_ulong", hwloc.}
## * \brief Remove index \p id from bitmap \p bitmap

proc hwloc_bitmap_clr*(bitmap: hwloc_bitmap_t; id: cuint): cint {.
    importc: "hwloc_bitmap_clr", hwloc.}
## * \brief Remove indexes from \p begin to \p end in bitmap \p bitmap.
##
##  If \p end is \c -1, the range is infinite.
##

proc hwloc_bitmap_clr_range*(bitmap: hwloc_bitmap_t; begin: cuint; `end`: cint): cint {.
    importc: "hwloc_bitmap_clr_range", hwloc.}
## * \brief Keep a single index among those set in bitmap \p bitmap
##
##  May be useful before binding so that the process does not
##  have a chance of migrating between multiple logical CPUs
##  in the original mask.
##  Instead of running the task on any PU inside the given CPU set,
##  the operating system scheduler will be forced to run it on a single
##  of these PUs.
##  It avoids a migration overhead and cache-line ping-pongs between PUs.
##
##  \note This function is NOT meant to distribute multiple processes
##  within a single CPU set. It always return the same single bit when
##  called multiple times on the same input set. hwloc_distrib() may
##  be used for generating CPU sets to distribute multiple tasks below
##  a single multi-PU object.
##
##  \note This function cannot be applied to an object set directly. It
##  should be applied to a copy (which may be obtained with hwloc_bitmap_dup()).
##

proc hwloc_bitmap_singlify*(bitmap: hwloc_bitmap_t): cint {.
    importc: "hwloc_bitmap_singlify", hwloc.}
##
##  Consulting bitmaps.
##
## * \brief Convert the beginning part of bitmap \p bitmap into unsigned long \p mask

proc hwloc_bitmap_to_ulong*(bitmap: hwloc_const_bitmap_t): culong {.
    importc: "hwloc_bitmap_to_ulong", hwloc.}
## * \brief Convert the \p i -th subset of bitmap \p bitmap into unsigned long mask

proc hwloc_bitmap_to_ith_ulong*(bitmap: hwloc_const_bitmap_t; i: cuint): culong {.
    importc: "hwloc_bitmap_to_ith_ulong", hwloc.}
## * \brief Convert the first \p nr subsets of bitmap \p bitmap into the array of \p nr unsigned long \p masks
##
##  \p nr may be determined earlier with hwloc_bitmap_nr_ulongs().
##
##  \return 0
##

proc hwloc_bitmap_to_ulongs*(bitmap: hwloc_const_bitmap_t; nr: cuint;
                            masks: ptr culong): cint {.
    importc: "hwloc_bitmap_to_ulongs", hwloc.}
## * \brief Return the number of unsigned longs required for storing bitmap \p bitmap entirely
##
##  This is the number of contiguous unsigned longs from the very first bit of the bitmap
##  (even if unset) up to the last set bit.
##  This is useful for knowing the \p nr parameter to pass to hwloc_bitmap_to_ulongs()
##  (or which calls to hwloc_bitmap_to_ith_ulong() are needed)
##  to entirely convert a bitmap into multiple unsigned longs.
##
##  When called on the output of hwloc_topology_get_topology_cpuset(),
##  the returned number is large enough for all cpusets of the topology.
##
##  \return -1 if \p bitmap is infinite.
##

proc hwloc_bitmap_nr_ulongs*(bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_nr_ulongs", hwloc.}
## * \brief Test whether index \p id is part of bitmap \p bitmap.
##
##  \return 1 if the bit at index \p id is set in bitmap \p bitmap, 0 otherwise.
##

proc hwloc_bitmap_isset*(bitmap: hwloc_const_bitmap_t; id: cuint): cint {.
    importc: "hwloc_bitmap_isset", hwloc.}
## * \brief Test whether bitmap \p bitmap is empty
##
##  \return 1 if bitmap is empty, 0 otherwise.
##

proc hwloc_bitmap_iszero*(bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_iszero", hwloc.}
## * \brief Test whether bitmap \p bitmap is completely full
##
##  \return 1 if bitmap is full, 0 otherwise.
##
##  \note A full bitmap is always infinitely set.
##

proc hwloc_bitmap_isfull*(bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_isfull", hwloc.}
## * \brief Compute the first index (least significant bit) in bitmap \p bitmap
##
##  \return -1 if no index is set in \p bitmap.
##

proc hwloc_bitmap_first*(bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_first", hwloc.}
## * \brief Compute the next index in bitmap \p bitmap which is after index \p prev
##
##  If \p prev is -1, the first index is returned.
##
##  \return -1 if no index with higher index is set in \p bitmap.
##

proc hwloc_bitmap_next*(bitmap: hwloc_const_bitmap_t; prev: cint): cint {.
    importc: "hwloc_bitmap_next", hwloc.}
## * \brief Compute the last index (most significant bit) in bitmap \p bitmap
##
##  \return -1 if no index is set in \p bitmap, or if \p bitmap is infinitely set.
##

proc hwloc_bitmap_last*(bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_last", hwloc.}
## * \brief Compute the "weight" of bitmap \p bitmap (i.e., number of
##  indexes that are in the bitmap).
##
##  \return the number of indexes that are in the bitmap.
##
##  \return -1 if \p bitmap is infinitely set.
##

proc hwloc_bitmap_weight*(bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_weight", hwloc.}
## * \brief Compute the first unset index (least significant bit) in bitmap \p bitmap
##
##  \return -1 if no index is unset in \p bitmap.
##

proc hwloc_bitmap_first_unset*(bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_first_unset", hwloc.}
## * \brief Compute the next unset index in bitmap \p bitmap which is after index \p prev
##
##  If \p prev is -1, the first unset index is returned.
##
##  \return -1 if no index with higher index is unset in \p bitmap.
##

proc hwloc_bitmap_next_unset*(bitmap: hwloc_const_bitmap_t; prev: cint): cint {.
    importc: "hwloc_bitmap_next_unset", hwloc.}
## * \brief Compute the last unset index (most significant bit) in bitmap \p bitmap
##
##  \return -1 if no index is unset in \p bitmap, or if \p bitmap is infinitely set.
##

proc hwloc_bitmap_last_unset*(bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_last_unset", hwloc.}
## * \brief Loop macro iterating on bitmap \p bitmap
##
##  The loop must start with hwloc_bitmap_foreach_begin() and end
##  with hwloc_bitmap_foreach_end() followed by a terminating ';'.
##
##  \p index is the loop variable; it should be an unsigned int.  The
##  first iteration will set \p index to the lowest index in the bitmap.
##  Successive iterations will iterate through, in order, all remaining
##  indexes set in the bitmap.  To be specific: each iteration will return a
##  value for \p index such that hwloc_bitmap_isset(bitmap, index) is true.
##
##  The assert prevents the loop from being infinite if the bitmap is infinitely set.
##
##  \hideinitializer
##
##  TODO: rewrite the macro in nim
##
##  Combining bitmaps.
##
## * \brief Or bitmaps \p bitmap1 and \p bitmap2 and store the result in bitmap \p res
##
##  \p res can be the same as \p bitmap1 or \p bitmap2
##

proc hwloc_bitmap_or*(res: hwloc_bitmap_t; bitmap1: hwloc_const_bitmap_t;
                     bitmap2: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_or", hwloc.}
## * \brief And bitmaps \p bitmap1 and \p bitmap2 and store the result in bitmap \p res
##
##  \p res can be the same as \p bitmap1 or \p bitmap2
##

proc hwloc_bitmap_and*(res: hwloc_bitmap_t; bitmap1: hwloc_const_bitmap_t;
                      bitmap2: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_and", hwloc.}
## * \brief And bitmap \p bitmap1 and the negation of \p bitmap2 and store the result in bitmap \p res
##
##  \p res can be the same as \p bitmap1 or \p bitmap2
##

proc hwloc_bitmap_andnot*(res: hwloc_bitmap_t; bitmap1: hwloc_const_bitmap_t;
                         bitmap2: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_andnot", hwloc.}
## * \brief Xor bitmaps \p bitmap1 and \p bitmap2 and store the result in bitmap \p res
##
##  \p res can be the same as \p bitmap1 or \p bitmap2
##

proc hwloc_bitmap_xor*(res: hwloc_bitmap_t; bitmap1: hwloc_const_bitmap_t;
                      bitmap2: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_xor", hwloc.}
## * \brief Negate bitmap \p bitmap and store the result in bitmap \p res
##
##  \p res can be the same as \p bitmap
##

proc hwloc_bitmap_not*(res: hwloc_bitmap_t; bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_not", hwloc.}
##
##  Comparing bitmaps.
##
## * \brief Test whether bitmaps \p bitmap1 and \p bitmap2 intersects.
##
##  \return 1 if bitmaps intersect, 0 otherwise.
##

proc hwloc_bitmap_intersects*(bitmap1: hwloc_const_bitmap_t;
                             bitmap2: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_intersects", hwloc.}
## * \brief Test whether bitmap \p sub_bitmap is part of bitmap \p super_bitmap.
##
##  \return 1 if \p sub_bitmap is included in \p super_bitmap, 0 otherwise.
##
##  \note The empty bitmap is considered included in any other bitmap.
##

proc hwloc_bitmap_isincluded*(sub_bitmap: hwloc_const_bitmap_t;
                             super_bitmap: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_isincluded", hwloc.}
## * \brief Test whether bitmap \p bitmap1 is equal to bitmap \p bitmap2.
##
##  \return 1 if bitmaps are equal, 0 otherwise.
##

proc hwloc_bitmap_isequal*(bitmap1: hwloc_const_bitmap_t;
                          bitmap2: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_isequal", hwloc.}
## * \brief Compare bitmaps \p bitmap1 and \p bitmap2 using their lowest index.
##
##  A bitmap is considered smaller if its least significant bit is smaller.
##  The empty bitmap is considered higher than anything (because its least significant bit does not exist).
##
##  \return -1 if \p bitmap1 is considered smaller than \p bitmap2.
##  \return 1 if \p bitmap1 is considered larger than \p bitmap2.
##
##  For instance comparing binary bitmaps 0011 and 0110 returns -1
##  (hence 0011 is considered smaller than 0110)
##  because least significant bit of 0011 (0001) is smaller than least significant bit of 0110 (0010).
##  Comparing 01001 and 00110 would also return -1 for the same reason.
##
##  \return 0 if bitmaps are considered equal, even if they are not strictly equal.
##  They just need to have the same least significant bit.
##  For instance, comparing binary bitmaps 0010 and 0110 returns 0 because they have the same least significant bit.
##

proc hwloc_bitmap_compare_first*(bitmap1: hwloc_const_bitmap_t;
                                bitmap2: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_compare_first", hwloc.}
## * \brief Compare bitmaps \p bitmap1 and \p bitmap2 in lexicographic order.
##
##  Lexicographic comparison of bitmaps, starting for their highest indexes.
##  Compare last indexes first, then second, etc.
##  The empty bitmap is considered lower than anything.
##
##  \return -1 if \p bitmap1 is considered smaller than \p bitmap2.
##  \return 1 if \p bitmap1 is considered larger than \p bitmap2.
##  \return 0 if bitmaps are equal (contrary to hwloc_bitmap_compare_first()).
##
##  For instance comparing binary bitmaps 0011 and 0110 returns -1
##  (hence 0011 is considered smaller than 0110).
##  Comparing 00101 and 01010 returns -1 too.
##
##  \note This is different from the non-existing hwloc_bitmap_compare_last()
##  which would only compare the highest index of each bitmap.
##

proc hwloc_bitmap_compare*(bitmap1: hwloc_const_bitmap_t;
                          bitmap2: hwloc_const_bitmap_t): cint {.
    importc: "hwloc_bitmap_compare", hwloc.}
## * @}

##
##  Copyright © 2009 CNRS
##  Copyright © 2009-2019 Inria.  All rights reserved.
##  Copyright © 2009-2012 Université Bordeaux
##  Copyright © 2009-2011 Cisco Systems, Inc.  All rights reserved.
##  See COPYING in top-level directory.
##
## =====================================================================
##                  PLEASE GO READ THE DOCUMENTATION!
##          ------------------------------------------------
##                $tarball_directory/doc/doxygen-doc/
##                                 or
##            http://www.open-mpi.org/projects/hwloc/doc/
## =====================================================================
##
##  FAIR WARNING: Do NOT expect to be able to figure out all the
##  subtleties of hwloc by simply reading function prototypes and
##  constant descrptions here in this file.
##
##  Hwloc has wonderful documentation in both PDF and HTML formats for
##  your reading pleasure.  The formal documentation explains a LOT of
##  hwloc-specific concepts, provides definitions, and discusses the
##  "big picture" for many of the things that you'll find here in this
##  header file.
##
##  The PDF/HTML documentation was generated via Doxygen; much of what
##  you'll see in there is also here in this file.  BUT THERE IS A LOT
##  THAT IS IN THE PDF/HTML THAT IS ***NOT*** IN hwloc.h!
##
##  There are entire paragraph-length descriptions, discussions, and
##  pretty prictures to explain subtle corner cases, provide concrete
##  examples, etc.
##
##  Please, go read the documentation.  :-)
##
##  Moreover there are several examples of hwloc use under doc/examples
##  in the source tree.
##
## =====================================================================
## * \file
##  \brief The hwloc API.
##
##  See hwloc/bitmap.h for bitmap specific macros.
##  See hwloc/helper.h for high-level topology traversal helpers.
##  See hwloc/inlines.h for the actual inline code of some functions below.
##  See hwloc/export.h for exporting topologies to XML or to synthetic descriptions.
##  See hwloc/distances.h for querying and modifying distances between objects.
##  See hwloc/diff.h for manipulating differences between similar topologies.
##

##
##  Symbol transforms
##

##
##  Bitmap definitions
##

## * \defgroup hwlocality_api_version API version
##  @{
##
## * \brief Indicate at build time which hwloc API version is being used.
##
##  This number is updated to (X<<16)+(Y<<8)+Z when a new release X.Y.Z
##  actually modifies the API.
##
##  Users may check for available features at build time using this number
##  (see \ref faq_version_api).
##
##  \note This should not be confused with HWLOC_VERSION, the library version.
##  Two stable releases of the same series usually have the same ::HWLOC_API_VERSION
##  even if their HWLOC_VERSION are different.
##

const
  HWLOC_API_VERSION* = 0x00020100

## * \brief Indicate at runtime which hwloc API version was used at build time.
##
##  Should be ::HWLOC_API_VERSION if running on the same version.
##

proc hwloc_get_api_version*(): cuint {.importc: "hwloc_get_api_version",
                                    hwloc.}
## * \brief Current component and plugin ABI version (see hwloc/plugins.h)

const
  HWLOC_COMPONENT_ABI* = 6

## * @}
## * \defgroup hwlocality_object_sets Object Sets (hwloc_cpuset_t and hwloc_nodeset_t)
##
##  Hwloc uses bitmaps to represent two distinct kinds of object sets:
##  CPU sets (::hwloc_cpuset_t) and NUMA node sets (::hwloc_nodeset_t).
##  These types are both typedefs to a common back end type
##  (::hwloc_bitmap_t), and therefore all the hwloc bitmap functions
##  are applicable to both ::hwloc_cpuset_t and ::hwloc_nodeset_t (see
##  \ref hwlocality_bitmap).
##
##  The rationale for having two different types is that even though
##  the actions one wants to perform on these types are the same (e.g.,
##  enable and disable individual items in the set/mask), they're used
##  in very different contexts: one for specifying which processors to
##  use and one for specifying which NUMA nodes to use.  Hence, the
##  name difference is really just to reflect the intent of where the
##  type is used.
##
##  @{
##
## * \brief A CPU set is a bitmap whose bits are set according to CPU
##  physical OS indexes.
##
##  It may be consulted and modified with the bitmap API as any
##  ::hwloc_bitmap_t (see hwloc/bitmap.h).
##
##  Each bit may be converted into a PU object using
##  hwloc_get_pu_obj_by_os_index().
##

type
  hwloc_cpuset_t* = hwloc_bitmap_t

## * \brief A non-modifiable ::hwloc_cpuset_t.

type
  hwloc_const_cpuset_t* = hwloc_const_bitmap_t

## * \brief A node set is a bitmap whose bits are set according to NUMA
##  memory node physical OS indexes.
##
##  It may be consulted and modified with the bitmap API as any
##  ::hwloc_bitmap_t (see hwloc/bitmap.h).
##  Each bit may be converted into a NUMA node object using
##  hwloc_get_numanode_obj_by_os_index().
##
##  When binding memory on a system without any NUMA node,
##  the single main memory bank is considered as NUMA node #0.
##
##  See also \ref hwlocality_helper_nodeset_convert.
##

type
  hwloc_nodeset_t* = hwloc_bitmap_t

## * \brief A non-modifiable ::hwloc_nodeset_t.
##

type
  hwloc_const_nodeset_t* = hwloc_const_bitmap_t

## * @}
## * \defgroup hwlocality_object_types Object Types
##  @{
##
## * \brief Type of topology object.
##
##  \note Do not rely on the ordering or completeness of the values as new ones
##  may be defined in the future!  If you need to compare types, use
##  hwloc_compare_types() instead.
##

type
  hwloc_obj_type_t* {.size: sizeof(cint).} = enum
    HWLOC_OBJ_MACHINE, ## *< \brief Machine.
                      ##  A set of processors and memory with cache
                      ##  coherency.
                      ##
                      ##  This type is always used for the root object of a topology,
                      ##  and never used anywhere else.
                      ##  Hence its parent is always \c NULL.
                      ##
    HWLOC_OBJ_PACKAGE, ## *< \brief Physical package.
                      ##  The physical package that usually gets inserted
                      ##  into a socket on the motherboard.
                      ##  A processor package usually contains multiple cores,
                      ##  and possibly some dies.
                      ##
    HWLOC_OBJ_CORE,           ## *< \brief Core.
                   ##  A computation unit (may be shared by several
                   ##  logical processors).
                   ##
    HWLOC_OBJ_PU, ## *< \brief Processing Unit, or (Logical) Processor.
                 ##  An execution unit (may share a core with some
                 ##  other logical processors, e.g. in the case of
                 ##  an SMT core).
                 ##
                 ##  This is the smallest object representing CPU resources,
                 ##  it cannot have any child except Misc objects.
                 ##
                 ##  Objects of this kind are always reported and can
                 ##  thus be used as fallback when others are not.
                 ##
    HWLOC_OBJ_L1CACHE,        ## *< \brief Level 1 Data (or Unified) Cache.
    HWLOC_OBJ_L2CACHE,        ## *< \brief Level 2 Data (or Unified) Cache.
    HWLOC_OBJ_L3CACHE,        ## *< \brief Level 3 Data (or Unified) Cache.
    HWLOC_OBJ_L4CACHE,        ## *< \brief Level 4 Data (or Unified) Cache.
    HWLOC_OBJ_L5CACHE,        ## *< \brief Level 5 Data (or Unified) Cache.
    HWLOC_OBJ_L1ICACHE,       ## *< \brief Level 1 instruction Cache (filtered out by default).
    HWLOC_OBJ_L2ICACHE,       ## *< \brief Level 2 instruction Cache (filtered out by default).
    HWLOC_OBJ_L3ICACHE,       ## *< \brief Level 3 instruction Cache (filtered out by default).
    HWLOC_OBJ_GROUP, ## *< \brief Group objects.
                    ##  Objects which do not fit in the above but are
                    ##  detected by hwloc and are useful to take into
                    ##  account for affinity. For instance, some operating systems
                    ##  expose their arbitrary processors aggregation this
                    ##  way.  And hwloc may insert such objects to group
                    ##  NUMA nodes according to their distances.
                    ##  See also \ref faq_groups.
                    ##
                    ##  These objects are removed when they do not bring
                    ##  any structure (see ::HWLOC_TYPE_FILTER_KEEP_STRUCTURE).
                    ##
    HWLOC_OBJ_NUMANODE, ## *< \brief NUMA node.
                       ##  An object that contains memory that is directly
                       ##  and byte-accessible to the host processors.
                       ##  It is usually close to some cores (the corresponding objects
                       ##  are descendants of the NUMA node object in the hwloc tree).
                       ##
                       ##  This is the smallest object representing Memory resources,
                       ##  it cannot have any child except Misc objects.
                       ##  However it may have Memory-side cache parents.
                       ##
                       ##  There is always at least one such object in the topology
                       ##  even if the machine is not NUMA.
                       ##
                       ##  Memory objects are not listed in the main children list,
                       ##  but rather in the dedicated Memory children list.
                       ##
                       ##  NUMA nodes have a special depth ::HWLOC_TYPE_DEPTH_NUMANODE
                       ##  instead of a normal depth just like other objects in the
                       ##  main tree.
                       ##
    HWLOC_OBJ_BRIDGE, ## *< \brief Bridge (filtered out by default).
                     ##  Any bridge that connects the host or an I/O bus,
                     ##  to another I/O bus.
                     ##  They are not added to the topology unless I/O discovery
                     ##  is enabled with hwloc_topology_set_flags().
                     ##  I/O objects are not listed in the main children list,
                     ##  but rather in the dedicated io children list.
                     ##  I/O objects have NULL CPU and node sets.
                     ##
    HWLOC_OBJ_PCI_DEVICE, ## *< \brief PCI device (filtered out by default).
                         ##  They are not added to the topology unless I/O discovery
                         ##  is enabled with hwloc_topology_set_flags().
                         ##  I/O objects are not listed in the main children list,
                         ##  but rather in the dedicated io children list.
                         ##  I/O objects have NULL CPU and node sets.
                         ##
    HWLOC_OBJ_OS_DEVICE, ## *< \brief Operating system device (filtered out by default).
                        ##  They are not added to the topology unless I/O discovery
                        ##  is enabled with hwloc_topology_set_flags().
                        ##  I/O objects are not listed in the main children list,
                        ##  but rather in the dedicated io children list.
                        ##  I/O objects have NULL CPU and node sets.
                        ##
    HWLOC_OBJ_MISC, ## *< \brief Miscellaneous objects (filtered out by default).
                   ##  Objects without particular meaning, that can e.g. be
                   ##  added by the application for its own use, or by hwloc
                   ##  for miscellaneous objects such as MemoryModule (DIMMs).
                   ##  These objects are not listed in the main children list,
                   ##  but rather in the dedicated misc children list.
                   ##  Misc objects may only have Misc objects as children,
                   ##  and those are in the dedicated misc children list as well.
                   ##  Misc objects have NULL CPU and node sets.
                   ##
    HWLOC_OBJ_MEMCACHE, ## *< \brief Memory-side cache (filtered out by default).
                       ##  A cache in front of a specific NUMA node.
                       ##
                       ##  This object always has at least one NUMA node as a memory child.
                       ##
                       ##  Memory objects are not listed in the main children list,
                       ##  but rather in the dedicated Memory children list.
                       ##
                       ##  Memory-side cache have a special depth ::HWLOC_TYPE_DEPTH_MEMCACHE
                       ##  instead of a normal depth just like other objects in the
                       ##  main tree.
                       ##
    HWLOC_OBJ_DIE, ## *< \brief Die within a physical package.
                  ##  A subpart of the physical package, that contains multiple cores.
                  ##  \hideinitializer
                  ##
    HWLOC_OBJ_TYPE_MAX        ## *< \private Sentinel value

const
  HWLOC_OBJ_TYPE_MIN* = HWLOC_OBJ_MACHINE


## * \brief Cache type.

type
  hwloc_obj_cache_type_t* {.size: sizeof(cint).} = enum
    HWLOC_OBJ_CACHE_UNIFIED,  ## *< \brief Unified cache.
    HWLOC_OBJ_CACHE_DATA,     ## *< \brief Data cache.
    HWLOC_OBJ_CACHE_INSTRUCTION ## *< \brief Instruction cache (filtered out by default).


## * \brief Type of one side (upstream or downstream) of an I/O bridge.

type
  hwloc_obj_bridge_type_t* {.size: sizeof(cint).} = enum
    HWLOC_OBJ_BRIDGE_HOST,    ## *< \brief Host-side of a bridge, only possible upstream.
    HWLOC_OBJ_BRIDGE_PCI      ## *< \brief PCI-side of a bridge.


## * \brief Type of a OS device.

type
  hwloc_obj_osdev_type_t* {.size: sizeof(cint).} = enum
    HWLOC_OBJ_OSDEV_BLOCK, ## *< \brief Operating system block device, or non-volatile memory device.
                          ##  For instance "sda" or "dax2.0" on Linux.
    HWLOC_OBJ_OSDEV_GPU,      ## *< \brief Operating system GPU device.
                        ##  For instance ":0.0" for a GL display,
                        ##  "card0" for a Linux DRM device.
    HWLOC_OBJ_OSDEV_NETWORK,  ## *< \brief Operating system network device.
                            ##  For instance the "eth0" interface on Linux.
    HWLOC_OBJ_OSDEV_OPENFABRICS, ## *< \brief Operating system openfabrics device.
                                ##  For instance the "mlx4_0" InfiniBand HCA,
                                ##  or "hfi1_0" Omni-Path interface on Linux.
    HWLOC_OBJ_OSDEV_DMA,      ## *< \brief Operating system dma engine device.
                        ##  For instance the "dma0chan0" DMA channel on Linux.
    HWLOC_OBJ_OSDEV_COPROC ## *< \brief Operating system co-processor device.
                          ##  For instance "mic0" for a Xeon Phi (MIC) on Linux,
                          ##  "opencl0d0" for a OpenCL device,
                          ##  "cuda0" for a CUDA device.


## * \brief Compare the depth of two object types
##
##  Types shouldn't be compared as they are, since newer ones may be added in
##  the future.  This function returns less than, equal to, or greater than zero
##  respectively if \p type1 objects usually include \p type2 objects, are the
##  same as \p type2 objects, or are included in \p type2 objects. If the types
##  can not be compared (because neither is usually contained in the other),
##  ::HWLOC_TYPE_UNORDERED is returned.  Object types containing CPUs can always
##  be compared (usually, a system contains machines which contain nodes which
##  contain packages which contain caches, which contain cores, which contain
##  processors).
##
##  \note ::HWLOC_OBJ_PU will always be the deepest,
##  while ::HWLOC_OBJ_MACHINE is always the highest.
##
##  \note This does not mean that the actual topology will respect that order:
##  e.g. as of today cores may also contain caches, and packages may also contain
##  nodes. This is thus just to be seen as a fallback comparison method.
##

proc hwloc_compare_types*(type1: hwloc_obj_type_t; type2: hwloc_obj_type_t): cint {.
    importc: "hwloc_compare_types", hwloc.}
type
  hwloc_compare_types_e* {.size: sizeof(cint).} = enum
    HWLOC_TYPE_UNORDERED = high(cint)


## * @}
## * \defgroup hwlocality_objects Object Structure and Attributes
##  @{
##

## * \brief Object type-specific Attributes

type
  hwloc_memory_page_type_s* {.bycopy.} = object
    size*: uint64      ## *< \brief Size of pages
    count*: uint64     ## *< \brief Number of pages of this size

  hwloc_numanode_attr_s* {.bycopy.} = object
    local_memory*: uint64 ## *< \brief Local memory (in bytes)
    page_types_len*: cuint ## *< \brief Size of array \p page_types
                         ## * \brief Array of local memory page types, \c NULL if no local memory and \p page_types is 0.
                         ##
                         ##  The array is sorted by increasing \p size fields.
                         ##  It contains \p page_types_len slots.
                         ##
    page_types*: ptr hwloc_memory_page_type_s

  hwloc_cache_attr_s* {.bycopy.} = object
    size*: uint64      ## *< \brief Size of cache in bytes
    depth*: cuint              ## *< \brief Depth of cache (e.g., L1, L2, ...etc.)
    linesize*: cuint           ## *< \brief Cache-line size in bytes. 0 if unknown
    associativity*: cint       ## *< \brief Ways of associativity,
                       ##   -1 if fully associative, 0 if unknown
    `type`*: hwloc_obj_cache_type_t ## *< \brief Cache type

  hwloc_group_attr_s* {.bycopy.} = object
    depth*: cuint              ## *< \brief Depth of group object.
                ##    It may change if intermediate Group objects are added.
    kind*: cuint               ## *< \brief Internally-used kind of group.
    subkind*: cuint            ## *< \brief Internally-used subkind to distinguish different levels of groups with same kind
    dont_merge*: uint8        ## *< \brief Flag preventing groups from being automatically merged with identical parent or children.

  hwloc_pcidev_attr_s* {.bycopy.} = object
    domain*: cushort
    bus*: uint8
    dev*: uint8
    `func`*: uint8
    class_id*: cushort
    vendor_id*: cushort
    device_id*: cushort
    subvendor_id*: cushort
    subdevice_id*: cushort
    revision*: uint8
    linkspeed*: cfloat         ##  in GB/s

  upstream_INNER_C_UNION* {.bycopy,union.} = object
    pci*: hwloc_pcidev_attr_s

  pci_INNER_C_STRUCT* {.bycopy.} = object
    domain*: cushort
    secondary_bus*: uint8
    subordinate_bus*: uint8

  downstream_INNER_C_UNION* {.bycopy,union.} = object
    pci*: pci_INNER_C_STRUCT

  hwloc_bridge_attr_s* {.bycopy.} = object
    upstream*: upstream_INNER_C_UNION
    upstream_type*: hwloc_obj_bridge_type_t
    downstream*: downstream_INNER_C_UNION
    downstream_type*: hwloc_obj_bridge_type_t
    depth*: cuint

  hwloc_osdev_attr_s* {.bycopy.} = object
    `type`*: hwloc_obj_osdev_type_t

  hwloc_obj_attr_u* {.bycopy,union.} = object
    numanode*: hwloc_numanode_attr_s ## * \brief NUMA node-specific Object Attributes
    ## * \brief Cache-specific Object Attributes
    cache*: hwloc_cache_attr_s ## * \brief Group-specific Object Attributes
    group*: hwloc_group_attr_s ## * \brief PCI Device specific Object Attributes
    pcidev*: hwloc_pcidev_attr_s ## * \brief Bridge specific Object Attribues
    bridge*: hwloc_bridge_attr_s ## * \brief OS Device specific Object Attributes
    osdev*: hwloc_osdev_attr_s


## * \brief Object info
##
##  \sa hwlocality_info_attr
##

type
  hwloc_info_s* {.bycopy.} = object
    name*: cstring             ## *< \brief Info name
    value*: cstring            ## *< \brief Info value


## * \brief Structure of a topology object
##
##  Applications must not modify any field except \p hwloc_obj.userdata.
##

const
  HWLOC_UNKNOWN_INDEX* = cast[cuint](cint(-1))

type
  hwloc_obj* {.bycopy.} = object
    `type`*: hwloc_obj_type_t  ##  physical information
    ## *< \brief Type of object
    subtype*: cstring          ## *< \brief Subtype string to better describe the type field.
    os_index*: cuint ## *< \brief OS-provided physical index number.
                   ##  It is not guaranteed unique across the entire machine,
                   ##  except for PUs and NUMA nodes.
                   ##  Set to HWLOC_UNKNOWN_INDEX if unknown or irrelevant for this object.
                   ##
    name*: cstring ## *< \brief Object-specific name if any.
                 ##  Mostly used for identifying OS devices and Misc objects where
                 ##  a name string is more useful than numerical indexes.
                 ##
    total_memory*: uint64 ## *< \brief Total memory (in bytes) in NUMA nodes below this object.
    attr*: ptr hwloc_obj_attr_u ## *< \brief Object type-specific Attributes,
                             ##  may be \c NULL if no attribute value was found
                             ##  global position
    depth*: cint ## *< \brief Vertical index in the hierarchy.
               ##
               ##  For normal objects, this is the depth of the horizontal level
               ##  that contains this object and its cousins of the same type.
               ##  If the topology is symmetric, this is equal to the parent depth
               ##  plus one, and also equal to the number of parent/child links
               ##  from the root object to here.
               ##
               ##  For special objects (NUMA nodes, I/O and Misc) that are not
               ##  in the main tree, this is a special negative value that
               ##  corresponds to their dedicated level,
               ##  see hwloc_get_type_depth() and ::hwloc_get_type_depth_e.
               ##  Those special values can be passed to hwloc functions such
               ##  hwloc_get_nbobjs_by_depth() as usual.
               ##
    logical_index*: cuint ## *< \brief Horizontal index in the whole list of similar objects,
                        ##  hence guaranteed unique across the entire machine.
                        ##  Could be a "cousin_rank" since it's the rank within the "cousin" list below
                        ##  Note that this index may change when restricting the topology
                        ##  or when inserting a group.
                        ##
                        ##  cousins are all objects of the same type (and depth) across the entire topology
    next_cousin*: ptr hwloc_obj ## *< \brief Next object of same type and depth
    prev_cousin*: ptr hwloc_obj ## *< \brief Previous object of same type and depth
                             ##  children of the same parent are siblings, even if they may have different type and depth
    parent*: ptr hwloc_obj      ## *< \brief Parent, \c NULL if root (Machine object)
    sibling_rank*: cuint       ## *< \brief Index in parent's \c children[] array. Or the index in parent's Memory, I/O or Misc children list.
    next_sibling*: ptr hwloc_obj ## *< \brief Next object below the same parent (inside the same list of children).
    prev_sibling*: ptr hwloc_obj ## *< \brief Previous object below the same parent (inside the same list of children).
                              ## * @name List and array of normal children below this object (except Memory, I/O and Misc children).
                              ## *@{
    arity*: cuint              ## *< \brief Number of normal children.
                ##  Memory, Misc and I/O children are not listed here
                ##  but rather in their dedicated children list.
                ##
    children*: ptr ptr hwloc_obj ## *< \brief Normal children, \c children[0 .. arity -1]
    first_child*: ptr hwloc_obj ## *< \brief First normal child
    last_child*: ptr hwloc_obj  ## *< \brief Last normal child
                            ## *@}
    symmetric_subtree*: cint ## *< \brief Set if the subtree of normal objects below this object is symmetric,
                           ##  which means all normal children and their children have identical subtrees.
                           ##
                           ##  Memory, I/O and Misc children are ignored.
                           ##
                           ##  If set in the topology root object, lstopo may export the topology
                           ##  as a synthetic string.
                           ##
                           ## * @name List of Memory children below this object.
                           ## *@{
    memory_arity*: cuint ## *< \brief Number of Memory children.
                       ##  These children are listed in \p memory_first_child.
                       ##
    memory_first_child*: ptr hwloc_obj ## *< \brief First Memory child.
                                    ##  NUMA nodes and Memory-side caches are listed here
                                    ##  (\p memory_arity and \p memory_first_child)
                                    ##  instead of in the normal children list.
                                    ##  See also hwloc_obj_type_is_memory().
                                    ##
                                    ##  A memory hierarchy starts from a normal CPU-side object
                                    ##  (e.g. Package) and ends with NUMA nodes as leaves.
                                    ##  There might exist some memory-side caches between them
                                    ##  in the middle of the memory subtree.
                                    ##
                                    ## *@}
                                    ## * @name List of I/O children below this object.
                                    ## *@{
    io_arity*: cuint           ## *< \brief Number of I/O children.
                   ##  These children are listed in \p io_first_child.
                   ##
    io_first_child*: ptr hwloc_obj ## *< \brief First I/O child.
                                ##  Bridges, PCI and OS devices are listed here (\p io_arity and \p io_first_child)
                                ##  instead of in the normal children list.
                                ##  See also hwloc_obj_type_is_io().
                                ##
                                ## *@}
                                ## * @name List of Misc children below this object.
                                ## *@{
    misc_arity*: cuint         ## *< \brief Number of Misc children.
                     ##  These children are listed in \p misc_first_child.
                     ##
    misc_first_child*: ptr hwloc_obj ## *< \brief First Misc child.
                                  ##  Misc objects are listed here (\p misc_arity and \p misc_first_child)
                                  ##  instead of in the normal children list.
                                  ##
                                  ## *@}
                                  ##  cpusets and nodesets
    cpuset*: hwloc_cpuset_t ## *< \brief CPUs covered by this object
                          ##
                          ##  This is the set of CPUs for which there are PU objects in the topology
                          ##  under this object, i.e. which are known to be physically contained in this
                          ##  object and known how (the children path between this object and the PU
                          ##  objects).
                          ##
                          ##  If the ::HWLOC_TOPOLOGY_FLAG_INCLUDE_DISALLOWED configuration flag is set,
                          ##  some of these CPUs may not be allowed for binding,
                          ##  see hwloc_topology_get_allowed_cpuset().
                          ##
                          ##  \note All objects have non-NULL CPU and node sets except Misc and I/O objects.
                          ##
                          ##  \note Its value must not be changed, hwloc_bitmap_dup() must be used instead.
                          ##
    complete_cpuset*: hwloc_cpuset_t ## *< \brief The complete CPU set of logical processors of this object,
                                   ##
                                   ##  This may include not only the same as the cpuset field, but also some CPUs for
                                   ##  which topology information is unknown or incomplete, some offlines CPUs, and
                                   ##  the CPUs that are ignored when the ::HWLOC_TOPOLOGY_FLAG_INCLUDE_DISALLOWED flag
                                   ##  is not set.
                                   ##  Thus no corresponding PU object may be found in the topology, because the
                                   ##  precise position is undefined. It is however known that it would be somewhere
                                   ##  under this object.
                                   ##
                                   ##  \note Its value must not be changed, hwloc_bitmap_dup() must be used instead.
                                   ##
    nodeset*: hwloc_nodeset_t ## *< \brief NUMA nodes covered by this object or containing this object
                            ##
                            ##  This is the set of NUMA nodes for which there are NUMA node objects in the
                            ##  topology under or above this object, i.e. which are known to be physically
                            ##  contained in this object or containing it and known how (the children path
                            ##  between this object and the NUMA node objects).
                            ##
                            ##  In the end, these nodes are those that are close to the current object.
                            ##
                            ##  If the ::HWLOC_TOPOLOGY_FLAG_INCLUDE_DISALLOWED configuration flag is set,
                            ##  some of these nodes may not be allowed for allocation,
                            ##  see hwloc_topology_get_allowed_nodeset().
                            ##
                            ##  If there are no NUMA nodes in the machine, all the memory is close to this
                            ##  object, so only the first bit may be set in \p nodeset.
                            ##
                            ##  \note All objects have non-NULL CPU and node sets except Misc and I/O objects.
                            ##
                            ##  \note Its value must not be changed, hwloc_bitmap_dup() must be used instead.
                            ##
    complete_nodeset*: hwloc_nodeset_t ## *< \brief The complete NUMA node set of this object,
                                     ##
                                     ##  This may include not only the same as the nodeset field, but also some NUMA
                                     ##  nodes for which topology information is unknown or incomplete, some offlines
                                     ##  nodes, and the nodes that are ignored when the ::HWLOC_TOPOLOGY_FLAG_INCLUDE_DISALLOWED
                                     ##  flag is not set.
                                     ##  Thus no corresponding NUMA node object may be found in the topology, because the
                                     ##  precise position is undefined. It is however known that it would be
                                     ##  somewhere under this object.
                                     ##
                                     ##  If there are no NUMA nodes in the machine, all the memory is close to this
                                     ##  object, so only the first bit is set in \p complete_nodeset.
                                     ##
                                     ##  \note Its value must not be changed, hwloc_bitmap_dup() must be used instead.
                                     ##
    infos*: ptr hwloc_info_s    ## *< \brief Array of stringified info type=name.
    infos_count*: cuint        ## *< \brief Size of infos array.
                      ##  misc
    userdata*: pointer ## *< \brief Application-given private data pointer,
                     ##  initialized to \c NULL, use it as you wish.
                     ##  See hwloc_topology_set_userdata_export_callback() in hwloc/export.h
                     ##  if you wish to export this field to XML.
    gp_index*: uint64 ## *< \brief Global persistent index.
                            ##  Generated by hwloc, unique across the topology (contrary to os_index)
                            ##  and persistent across topology changes (contrary to logical_index).
                            ##  Mostly used internally, but could also be used by application to identify objects.
                            ##


## *
##  \brief Convenience typedef; a pointer to a struct hwloc_obj.
##

type
  hwloc_obj_t* = ptr hwloc_obj

## * @}
## * \defgroup hwlocality_creation Topology Creation and Destruction
##  @{
##

type
  hwloc_topology* {.bycopy.} = object


## * \brief Topology context
##
##  To be initialized with hwloc_topology_init() and built with hwloc_topology_load().
##

type
  hwloc_topology_t* = ptr hwloc_topology

## * \brief Allocate a topology context.
##
##  \param[out] topologyp is assigned a pointer to the new allocated context.
##
##  \return 0 on success, -1 on error.
##

proc hwloc_topology_init*(topologyp: ptr hwloc_topology_t): cint {.
    importc: "hwloc_topology_init", hwloc.}
## * \brief Build the actual topology
##
##  Build the actual topology once initialized with hwloc_topology_init() and
##  tuned with \ref hwlocality_configuration and \ref hwlocality_setsource routines.
##  No other routine may be called earlier using this topology context.
##
##  \param topology is the topology to be loaded with objects.
##
##  \return 0 on success, -1 on error.
##
##  \note On failure, the topology is reinitialized. It should be either
##  destroyed with hwloc_topology_destroy() or configured and loaded again.
##
##  \note This function may be called only once per topology.
##
##  \note The binding of the current thread or process may temporarily change
##  during this call but it will be restored before it returns.
##
##  \sa hwlocality_configuration and hwlocality_setsource
##

proc hwloc_topology_load*(topology: hwloc_topology_t): cint {.
    importc: "hwloc_topology_load", hwloc.}
## * \brief Terminate and free a topology context
##
##  \param topology is the topology to be freed
##

proc hwloc_topology_destroy*(topology: hwloc_topology_t) {.
    importc: "hwloc_topology_destroy", hwloc.}
## * \brief Duplicate a topology.
##
##  The entire topology structure as well as its objects
##  are duplicated into a new one.
##
##  This is useful for keeping a backup while modifying a topology.
##
##  \note Object userdata is not duplicated since hwloc does not know what it point to.
##  The objects of both old and new topologies will point to the same userdata.
##

proc hwloc_topology_dup*(newtopology: ptr hwloc_topology_t;
                        oldtopology: hwloc_topology_t): cint {.
    importc: "hwloc_topology_dup", hwloc.}
## * \brief Verify that the topology is compatible with the current hwloc library.
##
##  This is useful when using the same topology structure (in memory)
##  in different libraries that may use different hwloc installations
##  (for instance if one library embeds a specific version of hwloc,
##  while another library uses a default system-wide hwloc installation).
##
##  If all libraries/programs use the same hwloc installation, this function
##  always returns success.
##
##  \return \c 0 on success.
##
##  \return \c -1 with \p errno set to \c EINVAL if incompatible.
##
##  \note If sharing between processes with hwloc_shmem_topology_write(),
##  the relevant check is already performed inside hwloc_shmem_topology_adopt().
##

proc hwloc_topology_abi_check*(topology: hwloc_topology_t): cint {.
    importc: "hwloc_topology_abi_check", hwloc.}
## * \brief Run internal checks on a topology structure
##
##  The program aborts if an inconsistency is detected in the given topology.
##
##  \param topology is the topology to be checked
##
##  \note This routine is only useful to developers.
##
##  \note The input topology should have been previously loaded with
##  hwloc_topology_load().
##

proc hwloc_topology_check*(topology: hwloc_topology_t) {.
    importc: "hwloc_topology_check", hwloc.}
## * @}
## * \defgroup hwlocality_levels Object levels, depths and types
##  @{
##
##  Be sure to see the figure in \ref termsanddefs that shows a
##  complete topology tree, including depths, child/sibling/cousin
##  relationships, and an example of an asymmetric topology where one
##  package has fewer caches than its peers.
##
## * \brief Get the depth of the hierarchical tree of objects.
##
##  This is the depth of ::HWLOC_OBJ_PU objects plus one.
##
##  \note NUMA nodes, I/O and Misc objects are ignored when computing
##  the depth of the tree (they are placed on special levels).
##

proc hwloc_topology_get_depth*(topology: hwloc_topology_t): cint {.
    importc: "hwloc_topology_get_depth", hwloc.}
## * \brief Returns the depth of objects of type \p type.
##
##  If no object of this type is present on the underlying architecture, or if
##  the OS doesn't provide this kind of information, the function returns
##  ::HWLOC_TYPE_DEPTH_UNKNOWN.
##
##  If type is absent but a similar type is acceptable, see also
##  hwloc_get_type_or_below_depth() and hwloc_get_type_or_above_depth().
##
##  If ::HWLOC_OBJ_GROUP is given, the function may return ::HWLOC_TYPE_DEPTH_MULTIPLE
##  if multiple levels of Groups exist.
##
##  If a NUMA node, I/O or Misc object type is given, the function returns a virtual
##  value because these objects are stored in special levels that are not CPU-related.
##  This virtual depth may be passed to other hwloc functions such as
##  hwloc_get_obj_by_depth() but it should not be considered as an actual
##  depth by the application. In particular, it should not be compared with
##  any other object depth or with the entire topology depth.
##  \sa hwloc_get_memory_parents_depth().
##
##  \sa hwloc_type_sscanf_as_depth() for returning the depth of objects
##  whose type is given as a string.
##

proc hwloc_get_type_depth*(topology: hwloc_topology_t; `type`: hwloc_obj_type_t): cint {.
    importc: "hwloc_get_type_depth", hwloc.}
type
  hwloc_get_type_depth_e* {.size: sizeof(cint).} = enum
    HWLOC_TYPE_DEPTH_MEMCACHE = -8, HWLOC_TYPE_DEPTH_MISC = -7, ## *< \brief Virtual depth for Misc object. \hideinitializer
    HWLOC_TYPE_DEPTH_OS_DEVICE = -6, ## *< \brief Virtual depth for software device object level. \hideinitializer
    HWLOC_TYPE_DEPTH_PCI_DEVICE = -5, ## *< \brief Virtual depth for PCI device object level. \hideinitializer
    HWLOC_TYPE_DEPTH_BRIDGE = -4, ## *< \brief Virtual depth for bridge object level. \hideinitializer
    HWLOC_TYPE_DEPTH_NUMANODE = -3, ## *< \brief Virtual depth for NUMA nodes. \hideinitializer
    HWLOC_TYPE_DEPTH_MULTIPLE = -2, ## *< \brief Objects of given type exist at different depth in the topology (only for Groups). \hideinitializer
    HWLOC_TYPE_DEPTH_UNKNOWN = -1 ## *< \brief No object of given type exists in the topology. \hideinitializer


## * \brief Return the depth of parents where memory objects are attached.
##
##  Memory objects have virtual negative depths because they are not part of
##  the main CPU-side hierarchy of objects. This depth should not be compared
##  with other level depths.
##
##  If all Memory objects are attached to Normal parents at the same depth,
##  this parent depth may be compared to other as usual, for instance
##  for knowing whether NUMA nodes is attached above or below Packages.
##
##  \return The depth of Normal parents of all memory children
##  if all these parents have the same depth. For instance the depth of
##  the Package level if all NUMA nodes are attached to Package objects.
##
##  \return ::HWLOC_TYPE_DEPTH_MULTIPLE if Normal parents of all
##  memory children do not have the same depth. For instance if some
##  NUMA nodes are attached to Packages while others are attached to
##  Groups.
##

proc hwloc_get_memory_parents_depth*(topology: hwloc_topology_t): cint {.
    importc: "hwloc_get_memory_parents_depth", hwloc.}
## * \brief Returns the depth of objects of type \p type or below
##
##  If no object of this type is present on the underlying architecture, the
##  function returns the depth of the first "present" object typically found
##  inside \p type.
##
##  This function is only meaningful for normal object types.
##  If a memory, I/O or Misc object type is given, the corresponding virtual
##  depth is always returned (see hwloc_get_type_depth()).
##
##  May return ::HWLOC_TYPE_DEPTH_MULTIPLE for ::HWLOC_OBJ_GROUP just like
##  hwloc_get_type_depth().
##

proc hwloc_get_type_or_below_depth*(topology: hwloc_topology_t;
                                   `type`: hwloc_obj_type_t): cint {.
    importc: "hwloc_get_type_or_below_depth", hwloc.}
## * \brief Returns the depth of objects of type \p type or above
##
##  If no object of this type is present on the underlying architecture, the
##  function returns the depth of the first "present" object typically
##  containing \p type.
##
##  This function is only meaningful for normal object types.
##  If a memory, I/O or Misc object type is given, the corresponding virtual
##  depth is always returned (see hwloc_get_type_depth()).
##
##  May return ::HWLOC_TYPE_DEPTH_MULTIPLE for ::HWLOC_OBJ_GROUP just like
##  hwloc_get_type_depth().
##

proc hwloc_get_type_or_above_depth*(topology: hwloc_topology_t;
                                   `type`: hwloc_obj_type_t): cint {.
    importc: "hwloc_get_type_or_above_depth", hwloc.}
## * \brief Returns the type of objects at depth \p depth.
##
##  \p depth should between 0 and hwloc_topology_get_depth()-1.
##
##  \return (hwloc_obj_type_t)-1 if depth \p depth does not exist.
##

proc hwloc_get_depth_type*(topology: hwloc_topology_t; depth: cint): hwloc_obj_type_t {.
    importc: "hwloc_get_depth_type", hwloc.}
## * \brief Returns the width of level at depth \p depth.
##

proc hwloc_get_nbobjs_by_depth*(topology: hwloc_topology_t; depth: cint): cuint {.
    importc: "hwloc_get_nbobjs_by_depth", hwloc.}
## * \brief Returns the width of level type \p type
##
##  If no object for that type exists, 0 is returned.
##  If there are several levels with objects of that type, -1 is returned.
##

proc hwloc_get_nbobjs_by_type*(topology: hwloc_topology_t; `type`: hwloc_obj_type_t): cint {.
    importc: "hwloc_get_nbobjs_by_type", hwloc.}
## * \brief Returns the top-object of the topology-tree.
##
##  Its type is ::HWLOC_OBJ_MACHINE.
##

proc hwloc_get_root_obj*(topology: hwloc_topology_t): hwloc_obj_t {.
    importc: "hwloc_get_root_obj", hwloc.}
## * \brief Returns the topology object at logical index \p idx from depth \p depth

proc hwloc_get_obj_by_depth*(topology: hwloc_topology_t; depth: cint; idx: cuint): hwloc_obj_t {.
    importc: "hwloc_get_obj_by_depth", hwloc.}
## * \brief Returns the topology object at logical index \p idx with type \p type
##
##  If no object for that type exists, \c NULL is returned.
##  If there are several levels with objects of that type (::HWLOC_OBJ_GROUP),
##  \c NULL is returned and the caller may fallback to hwloc_get_obj_by_depth().
##

proc hwloc_get_obj_by_type*(topology: hwloc_topology_t; `type`: hwloc_obj_type_t;
                           idx: cuint): hwloc_obj_t {.
    importc: "hwloc_get_obj_by_type", hwloc.}
## * \brief Returns the next object at depth \p depth.
##
##  If \p prev is \c NULL, return the first object at depth \p depth.
##

proc hwloc_get_next_obj_by_depth*(topology: hwloc_topology_t; depth: cint;
                                 prev: hwloc_obj_t): hwloc_obj_t {.
    importc: "hwloc_get_next_obj_by_depth", hwloc.}
## * \brief Returns the next object of type \p type.
##
##  If \p prev is \c NULL, return the first object at type \p type.  If
##  there are multiple or no depth for given type, return \c NULL and
##  let the caller fallback to hwloc_get_next_obj_by_depth().
##

proc hwloc_get_next_obj_by_type*(topology: hwloc_topology_t;
                                `type`: hwloc_obj_type_t; prev: hwloc_obj_t): hwloc_obj_t {.
    importc: "hwloc_get_next_obj_by_type", hwloc.}
## * @}
## * \defgroup hwlocality_object_strings Converting between Object Types and Attributes, and Strings
##  @{
##
## * \brief Return a constant stringified object type.
##
##  This function is the basic way to convert a generic type into a string.
##  The output string may be parsed back by hwloc_type_sscanf().
##
##  hwloc_obj_type_snprintf() may return a more precise output for a specific
##  object, but it requires the caller to provide the output buffer.
##

proc hwloc_obj_type_string*(`type`: hwloc_obj_type_t): cstring {.
    importc: "hwloc_obj_type_string", hwloc.}
## * \brief Stringify the type of a given topology object into a human-readable form.
##
##  Contrary to hwloc_obj_type_string(), this function includes object-specific
##  attributes (such as the Group depth, the Bridge type, or OS device type)
##  in the output, and it requires the caller to provide the output buffer.
##
##  The output is guaranteed to be the same for all objects of a same topology level.
##
##  If \p verbose is 1, longer type names are used, e.g. L1Cache instead of L1.
##
##  The output string may be parsed back by hwloc_type_sscanf().
##
##  If \p size is 0, \p string may safely be \c NULL.
##
##  \return the number of character that were actually written if not truncating,
##  or that would have been written (not including the ending \\0).
##

proc hwloc_obj_type_snprintf*(string: cstring; size: csize_t; obj: hwloc_obj_t;
                             verbose: cint): cint {.
    importc: "hwloc_obj_type_snprintf", hwloc.}
## * \brief Stringify the attributes of a given topology object into a human-readable form.
##
##  Attribute values are separated by \p separator.
##
##  Only the major attributes are printed in non-verbose mode.
##
##  If \p size is 0, \p string may safely be \c NULL.
##
##  \return the number of character that were actually written if not truncating,
##  or that would have been written (not including the ending \\0).
##

proc hwloc_obj_attr_snprintf*(string: cstring; size: csize_t; obj: hwloc_obj_t;
                             separator: cstring; verbose: cint): cint {.
    importc: "hwloc_obj_attr_snprintf", hwloc.}
## * \brief Return an object type and attributes from a type string.
##
##  Convert strings such as "Package" or "L1iCache" into the corresponding types.
##  Matching is case-insensitive, and only the first letters are actually
##  required to match.
##
##  The matched object type is set in \p typep (which cannot be \c NULL).
##
##  Type-specific attributes, for instance Cache type, Cache depth, Group depth,
##  Bridge type or OS Device type may be returned in \p attrp.
##  Attributes that are not specified in the string (for instance "Group"
##  without a depth, or "L2Cache" without a cache type) are set to -1.
##
##  \p attrp is only filled if not \c NULL and if its size specified in \p attrsize
##  is large enough. It should be at least as large as union hwloc_obj_attr_u.
##
##  \return 0 if a type was correctly identified, otherwise -1.
##
##  \note This function is guaranteed to match any string returned by
##  hwloc_obj_type_string() or hwloc_obj_type_snprintf().
##
##  \note This is an extended version of the now deprecated hwloc_obj_type_sscanf().
##

proc hwloc_type_sscanf*(string: cstring; typep: ptr hwloc_obj_type_t;
                       attrp: ptr hwloc_obj_attr_u; attrsize: csize_t): cint {.
    importc: "hwloc_type_sscanf", hwloc.}
## * \brief Return an object type and its level depth from a type string.
##
##  Convert strings such as "Package" or "L1iCache" into the corresponding types
##  and return in \p depthp the depth of the corresponding level in the
##  topology \p topology.
##
##  If no object of this type is present on the underlying architecture,
##  ::HWLOC_TYPE_DEPTH_UNKNOWN is returned.
##
##  If multiple such levels exist (for instance if giving Group without any depth),
##  the function may return ::HWLOC_TYPE_DEPTH_MULTIPLE instead.
##
##  The matched object type is set in \p typep if \p typep is non \c NULL.
##
##  \note This function is similar to hwloc_type_sscanf() followed
##  by hwloc_get_type_depth() but it also automatically disambiguates
##  multiple group levels etc.
##
##  \note This function is guaranteed to match any string returned by
##  hwloc_obj_type_string() or hwloc_obj_type_snprintf().
##

proc hwloc_type_sscanf_as_depth*(string: cstring; typep: ptr hwloc_obj_type_t;
                                topology: hwloc_topology_t; depthp: ptr cint): cint {.
    importc: "hwloc_type_sscanf_as_depth", hwloc.}
## * @}
## * \defgroup hwlocality_info_attr Consulting and Adding Key-Value Info Attributes
##
##  @{
##
## * \brief Search the given key name in object infos and return the corresponding value.
##
##  If multiple keys match the given name, only the first one is returned.
##
##  \return \c NULL if no such key exists.
##

proc hwloc_obj_get_info_by_name*(obj: hwloc_obj_t; name: cstring): cstring {.
    importc: "hwloc_obj_get_info_by_name", hwloc.}
## * \brief Add the given info name and value pair to the given object.
##
##  The info is appended to the existing info array even if another key
##  with the same name already exists.
##
##  The input strings are copied before being added in the object infos.
##
##  \return \c 0 on success, \c -1 on error.
##
##  \note This function may be used to enforce object colors in the lstopo
##  graphical output by using "lstopoStyle" as a name and "Background=#rrggbb"
##  as a value. See CUSTOM COLORS in the lstopo(1) manpage for details.
##
##  \note If \p value contains some non-printable characters, they will
##  be dropped when exporting to XML, see hwloc_topology_export_xml() in hwloc/export.h.
##

proc hwloc_obj_add_info*(obj: hwloc_obj_t; name: cstring; value: cstring): cint {.
    importc: "hwloc_obj_add_info", hwloc.}
## * @}
## * \defgroup hwlocality_cpubinding CPU binding
##
##  Some operating systems only support binding threads or processes to a single PU.
##  Others allow binding to larger sets such as entire Cores or Packages or
##  even random sets of invididual PUs. In such operating system, the scheduler
##  is free to run the task on one of these PU, then migrate it to another PU, etc.
##  It is often useful to call hwloc_bitmap_singlify() on the target CPU set before
##  passing it to the binding function to avoid these expensive migrations.
##  See the documentation of hwloc_bitmap_singlify() for details.
##
##  Some operating systems do not provide all hwloc-supported
##  mechanisms to bind processes, threads, etc.
##  hwloc_topology_get_support() may be used to query about the actual CPU
##  binding support in the currently used operating system.
##
##  When the requested binding operation is not available and the
##  ::HWLOC_CPUBIND_STRICT flag was passed, the function returns -1.
##  \p errno is set to \c ENOSYS when it is not possible to bind the requested kind of object
##  processes/threads. errno is set to \c EXDEV when the requested cpuset
##  can not be enforced (e.g. some systems only allow one CPU, and some
##  other systems only allow one NUMA node).
##
##  If ::HWLOC_CPUBIND_STRICT was not passed, the function may fail as well,
##  or the operating system may use a slightly different operation
##  (with side-effects, smaller binding set, etc.)
##  when the requested operation is not exactly supported.
##
##  The most portable version that should be preferred over the others,
##  whenever possible, is the following one which just binds the current program,
##  assuming it is single-threaded:
##
##  \code
##  hwloc_set_cpubind(topology, set, 0),
##  \endcode
##
##  If the program may be multithreaded, the following one should be preferred
##  to only bind the current thread:
##
##  \code
##  hwloc_set_cpubind(topology, set, HWLOC_CPUBIND_THREAD),
##  \endcode
##
##  \sa Some example codes are available under doc/examples/ in the source tree.
##
##  \note To unbind, just call the binding function with either a full cpuset or
##  a cpuset equal to the system cpuset.
##
##  \note On some operating systems, CPU binding may have effects on memory binding, see
##  ::HWLOC_CPUBIND_NOMEMBIND
##
##  \note Running lstopo \--top or hwloc-ps can be a very convenient tool to check
##  how binding actually happened.
##  @{
##
## * \brief Process/Thread binding flags.
##
##  These bit flags can be used to refine the binding policy.
##
##  The default (0) is to bind the current process, assumed to be
##  single-threaded, in a non-strict way.  This is the most portable
##  way to bind as all operating systems usually provide it.
##
##  \note Not all systems support all kinds of binding.  See the
##  "Detailed Description" section of \ref hwlocality_cpubinding for a
##  description of errors that can occur.
##

type ## * \brief Bind all threads of the current (possibly) multithreaded process.
    ##  \hideinitializer
  hwloc_cpubind_flags_t* {.size: sizeof(cint).} = enum
    HWLOC_CPUBIND_PROCESS = (1 shl 0), ## * \brief Bind current thread of current process.
                                  ##  \hideinitializer
    HWLOC_CPUBIND_THREAD = (1 shl 1), ## * \brief Request for strict binding from the OS.
                                 ##
                                 ##  By default, when the designated CPUs are all busy while other
                                 ##  CPUs are idle, operating systems may execute the thread/process
                                 ##  on those other CPUs instead of the designated CPUs, to let them
                                 ##  progress anyway.  Strict binding means that the thread/process
                                 ##  will _never_ execute on other cpus than the designated CPUs, even
                                 ##  when those are busy with other tasks and other CPUs are idle.
                                 ##
                                 ##  \note Depending on the operating system, strict binding may not
                                 ##  be possible (e.g., the OS does not implement it) or not allowed
                                 ##  (e.g., for an administrative reasons), and the function will fail
                                 ##  in that case.
                                 ##
                                 ##  When retrieving the binding of a process, this flag checks
                                 ##  whether all its threads  actually have the same binding. If the
                                 ##  flag is not given, the binding of each thread will be
                                 ##  accumulated.
                                 ##
                                 ##  \note This flag is meaningless when retrieving the binding of a
                                 ##  thread.
                                 ##  \hideinitializer
                                 ##
    HWLOC_CPUBIND_STRICT = (1 shl 2), ## * \brief Avoid any effect on memory binding
                                 ##
                                 ##  On some operating systems, some CPU binding function would also
                                 ##  bind the memory on the corresponding NUMA node.  It is often not
                                 ##  a problem for the application, but if it is, setting this flag
                                 ##  will make hwloc avoid using OS functions that would also bind
                                 ##  memory.  This will however reduce the support of CPU bindings,
                                 ##  i.e. potentially return -1 with errno set to ENOSYS in some
                                 ##  cases.
                                 ##
                                 ##  This flag is only meaningful when used with functions that set
                                 ##  the CPU binding.  It is ignored when used with functions that get
                                 ##  CPU binding information.
                                 ##  \hideinitializer
                                 ##
    HWLOC_CPUBIND_NOMEMBIND = (1 shl 3)


from posix import Pid, Pthread
type
  hwloc_pid_t* = Pid
  hwloc_thread_t* = Pthread

## * \brief Bind current process or thread on cpus given in physical bitmap \p set.
##
##  \return -1 with errno set to ENOSYS if the action is not supported
##  \return -1 with errno set to EXDEV if the binding cannot be enforced
##

proc hwloc_set_cpubind*(topology: hwloc_topology_t; set: hwloc_const_cpuset_t;
                       flags: cint): cint {.importc: "hwloc_set_cpubind",
    hwloc.}
## * \brief Get current process or thread binding.
##
##  Writes into \p set the physical cpuset which the process or thread (according to \e
##  flags) was last bound to.
##

proc hwloc_get_cpubind*(topology: hwloc_topology_t; set: hwloc_cpuset_t; flags: cint): cint {.
    importc: "hwloc_get_cpubind", hwloc.}
## * \brief Bind a process \p pid on cpus given in physical bitmap \p set.
##
##  \note \p hwloc_pid_t is \p pid_t on Unix platforms,
##  and \p HANDLE on native Windows platforms.
##
##  \note As a special case on Linux, if a tid (thread ID) is supplied
##  instead of a pid (process ID) and ::HWLOC_CPUBIND_THREAD is passed in flags,
##  the binding is applied to that specific thread.
##
##  \note On non-Linux systems, ::HWLOC_CPUBIND_THREAD can not be used in \p flags.
##

proc hwloc_set_proc_cpubind*(topology: hwloc_topology_t; pid: hwloc_pid_t;
                            set: hwloc_const_cpuset_t; flags: cint): cint {.
    importc: "hwloc_set_proc_cpubind", hwloc.}
## * \brief Get the current physical binding of process \p pid.
##
##  \note \p hwloc_pid_t is \p pid_t on Unix platforms,
##  and \p HANDLE on native Windows platforms.
##
##  \note As a special case on Linux, if a tid (thread ID) is supplied
##  instead of a pid (process ID) and HWLOC_CPUBIND_THREAD is passed in flags,
##  the binding for that specific thread is returned.
##
##  \note On non-Linux systems, HWLOC_CPUBIND_THREAD can not be used in \p flags.
##

proc hwloc_get_proc_cpubind*(topology: hwloc_topology_t; pid: hwloc_pid_t;
                            set: hwloc_cpuset_t; flags: cint): cint {.
    importc: "hwloc_get_proc_cpubind", hwloc.}
## * \brief Bind a thread \p thread on cpus given in physical bitmap \p set.
##
##  \note \p hwloc_thread_t is \p pthread_t on Unix platforms,
##  and \p HANDLE on native Windows platforms.
##
##  \note ::HWLOC_CPUBIND_PROCESS can not be used in \p flags.
##
proc hwloc_set_thread_cpubind*(topology: hwloc_topology_t;
                              thread: hwloc_thread_t; set: hwloc_const_cpuset_t;
                              flags: cint): cint {.
    importc: "hwloc_set_thread_cpubind", hwloc.}
## * \brief Get the current physical binding of thread \p tid.
##
##  \note \p hwloc_thread_t is \p pthread_t on Unix platforms,
##  and \p HANDLE on native Windows platforms.
##
##  \note ::HWLOC_CPUBIND_PROCESS can not be used in \p flags.
##
proc hwloc_get_thread_cpubind*(topology: hwloc_topology_t;
                              thread: hwloc_thread_t; set: hwloc_cpuset_t;
                              flags: cint): cint {.
    importc: "hwloc_get_thread_cpubind", hwloc.}
## * \brief Get the last physical CPU where the current process or thread ran.
##
##  The operating system may move some tasks from one processor
##  to another at any time according to their binding,
##  so this function may return something that is already
##  outdated.
##
##  \p flags can include either ::HWLOC_CPUBIND_PROCESS or ::HWLOC_CPUBIND_THREAD to
##  specify whether the query should be for the whole process (union of all CPUs
##  on which all threads are running), or only the current thread. If the
##  process is single-threaded, flags can be set to zero to let hwloc use
##  whichever method is available on the underlying OS.
##

proc hwloc_get_last_cpu_location*(topology: hwloc_topology_t; set: hwloc_cpuset_t;
                                 flags: cint): cint {.
    importc: "hwloc_get_last_cpu_location", hwloc.}
## * \brief Get the last physical CPU where a process ran.
##
##  The operating system may move some tasks from one processor
##  to another at any time according to their binding,
##  so this function may return something that is already
##  outdated.
##
##  \note \p hwloc_pid_t is \p pid_t on Unix platforms,
##  and \p HANDLE on native Windows platforms.
##
##  \note As a special case on Linux, if a tid (thread ID) is supplied
##  instead of a pid (process ID) and ::HWLOC_CPUBIND_THREAD is passed in flags,
##  the last CPU location of that specific thread is returned.
##
##  \note On non-Linux systems, ::HWLOC_CPUBIND_THREAD can not be used in \p flags.
##

proc hwloc_get_proc_last_cpu_location*(topology: hwloc_topology_t;
                                      pid: hwloc_pid_t; set: hwloc_cpuset_t;
                                      flags: cint): cint {.
    importc: "hwloc_get_proc_last_cpu_location", hwloc.}
## * @}
## * \defgroup hwlocality_membinding Memory binding
##
##  Memory binding can be done three ways:
##
##  - explicit memory allocation thanks to hwloc_alloc_membind() and friends:
##    the binding will have effect on the memory allocated by these functions.
##  - implicit memory binding through binding policy: hwloc_set_membind() and
##    friends only define the current policy of the process, which will be
##    applied to the subsequent calls to malloc() and friends.
##  - migration of existing memory ranges, thanks to hwloc_set_area_membind()
##    and friends, which move already-allocated data.
##
##  Not all operating systems support all three ways.
##  hwloc_topology_get_support() may be used to query about the actual memory
##  binding support in the currently used operating system.
##
##  When the requested binding operation is not available and the
##  ::HWLOC_MEMBIND_STRICT flag was passed, the function returns -1.
##  \p errno will be set to \c ENOSYS when the system does support
##  the specified action or policy
##  (e.g., some systems only allow binding memory on a per-thread
##  basis, whereas other systems only allow binding memory for all
##  threads in a process).
##  \p errno will be set to EXDEV when the requested set can not be enforced
##  (e.g., some systems only allow binding memory to a single NUMA node).
##
##  If ::HWLOC_MEMBIND_STRICT was not passed, the function may fail as well,
##  or the operating system may use a slightly different operation
##  (with side-effects, smaller binding set, etc.)
##  when the requested operation is not exactly supported.
##
##  The most portable form that should be preferred over the others
##  whenever possible is as follows.
##  It allocates some memory hopefully bound to the specified set.
##  To do so, hwloc will possibly have to change the current memory
##  binding policy in order to actually get the memory bound, if the OS
##  does not provide any other way to simply allocate bound memory
##  without changing the policy for all allocations. That is the
##  difference with hwloc_alloc_membind(), which will never change the
##  current memory binding policy.
##
##  \code
##  hwloc_alloc_membind_policy(topology, size, set,
##                             HWLOC_MEMBIND_BIND, 0);
##  \endcode
##
##  Each hwloc memory binding function takes a bitmap argument that
##  is a CPU set by default, or a NUMA memory node set if the flag
##  ::HWLOC_MEMBIND_BYNODESET is specified.
##  See \ref hwlocality_object_sets and \ref hwlocality_bitmap for a
##  discussion of CPU sets and NUMA memory node sets.
##  It is also possible to convert between CPU set and node set using
##  hwloc_cpuset_to_nodeset() or hwloc_cpuset_from_nodeset().
##
##  Memory binding by CPU set cannot work for CPU-less NUMA memory nodes.
##  Binding by nodeset should therefore be preferred whenever possible.
##
##  \sa Some example codes are available under doc/examples/ in the source tree.
##
##  \note On some operating systems, memory binding affects the CPU
##  binding; see ::HWLOC_MEMBIND_NOCPUBIND
##  @{
##
## * \brief Memory binding policy.
##
##  These constants can be used to choose the binding policy.  Only one policy can
##  be used at a time (i.e., the values cannot be OR'ed together).
##
##  Not all systems support all kinds of binding.
##  hwloc_topology_get_support() may be used to query about the actual memory
##  binding policy support in the currently used operating system.
##  See the "Detailed Description" section of \ref hwlocality_membinding
##  for a description of errors that can occur.
##

type ## * \brief Reset the memory allocation policy to the system default.
    ##  Depending on the operating system, this may correspond to
    ##  ::HWLOC_MEMBIND_FIRSTTOUCH (Linux),
    ##  or ::HWLOC_MEMBIND_BIND (AIX, HP-UX, Solaris, Windows).
    ##  This policy is never returned by get membind functions.
    ##  The nodeset argument is ignored.
    ##  \hideinitializer
  hwloc_membind_policy_t* {.size: sizeof(cint).} = enum
    HWLOC_MEMBIND_MIXED = -1, HWLOC_MEMBIND_DEFAULT = 0, ## * \brief Allocate each memory page individually on the local NUMA
                                                   ##  node of the thread that touches it.
                                                   ##
                                                   ##  The given nodeset should usually be
                                                   ## hwloc_topology_get_topology_nodeset()
                                                   ##  so that the touching thread may run and allocate on any node in the system.
                                                   ##
                                                   ##  On AIX, if the nodeset is smaller, pages are allocated locally (if the local
                                                   ##  node is in the nodeset) or from a random non-local node (otherwise).
                                                   ##  \hideinitializer
    HWLOC_MEMBIND_FIRSTTOUCH = 1, ## * \brief Allocate memory on the specified nodes.
                               ##  \hideinitializer
    HWLOC_MEMBIND_BIND = 2, ## * \brief Allocate memory on the given nodes in an interleaved
                         ##  / round-robin manner.  The precise layout of the memory across
                         ##  multiple NUMA nodes is OS/system specific. Interleaving can be
                         ##  useful when threads distributed across the specified NUMA nodes
                         ##  will all be accessing the whole memory range concurrently, since
                         ##  the interleave will then balance the memory references.
                         ##  \hideinitializer
    HWLOC_MEMBIND_INTERLEAVE = 3, ## * \brief For each page bound with this policy, by next time
                               ##  it is touched (and next time only), it is moved from its current
                               ##  location to the local NUMA node of the thread where the memory
                               ##  reference occurred (if it needs to be moved at all).
                               ##  \hideinitializer
    HWLOC_MEMBIND_NEXTTOUCH = 4 ## * \brief Returned by get_membind() functions when multiple
                             ##  threads or parts of a memory area have differing memory binding
                             ##  policies.
                             ##  Also returned when binding is unknown because binding hooks are empty
                             ##  when the topology is loaded from XML without HWLOC_THISSYSTEM=1, etc.
                             ##  \hideinitializer


## * \brief Memory binding flags.
##
##  These flags can be used to refine the binding policy.
##  All flags can be logically OR'ed together with the exception of
##  ::HWLOC_MEMBIND_PROCESS and ::HWLOC_MEMBIND_THREAD;
##  these two flags are mutually exclusive.
##
##  Not all systems support all kinds of binding.
##  hwloc_topology_get_support() may be used to query about the actual memory
##  binding support in the currently used operating system.
##  See the "Detailed Description" section of \ref hwlocality_membinding
##  for a description of errors that can occur.
##

type ## * \brief Set policy for all threads of the specified (possibly
    ##  multithreaded) process.  This flag is mutually exclusive with
    ##  ::HWLOC_MEMBIND_THREAD.
    ##  \hideinitializer
  hwloc_membind_flags_t* {.size: sizeof(cint).} = enum
    HWLOC_MEMBIND_PROCESS = (1 shl 0), ## * \brief Set policy for a specific thread of the current process.
                                  ##  This flag is mutually exclusive with ::HWLOC_MEMBIND_PROCESS.
                                  ##  \hideinitializer
    HWLOC_MEMBIND_THREAD = (1 shl 1), ## * Request strict binding from the OS.  The function will fail if
                                 ##  the binding can not be guaranteed / completely enforced.
                                 ##
                                 ##  This flag has slightly different meanings depending on which
                                 ##  function it is used with.
                                 ##  \hideinitializer
    HWLOC_MEMBIND_STRICT = (1 shl 2), ## * \brief Migrate existing allocated memory.  If the memory cannot
                                 ##  be migrated and the ::HWLOC_MEMBIND_STRICT flag is passed, an error
                                 ##  will be returned.
                                 ##  \hideinitializer
    HWLOC_MEMBIND_MIGRATE = (1 shl 3), ## * \brief Avoid any effect on CPU binding.
                                  ##
                                  ##  On some operating systems, some underlying memory binding
                                  ##  functions also bind the application to the corresponding CPU(s).
                                  ##  Using this flag will cause hwloc to avoid using OS functions that
                                  ##  could potentially affect CPU bindings.  Note, however, that using
                                  ##  NOCPUBIND may reduce hwloc's overall memory binding
                                  ##  support. Specifically: some of hwloc's memory binding functions
                                  ##  may fail with errno set to ENOSYS when used with NOCPUBIND.
                                  ##  \hideinitializer
                                  ##
    HWLOC_MEMBIND_NOCPUBIND = (1 shl 4), ## * \brief Consider the bitmap argument as a nodeset.
                                    ##
                                    ##  The bitmap argument is considered a nodeset if this flag is given,
                                    ##  or a cpuset otherwise by default.
                                    ##
                                    ##  Memory binding by CPU set cannot work for CPU-less NUMA memory nodes.
                                    ##  Binding by nodeset should therefore be preferred whenever possible.
                                    ##  \hideinitializer
                                    ##
    HWLOC_MEMBIND_BYNODESET = (1 shl 5)


## * \brief Set the default memory binding policy of the current
##  process or thread to prefer the NUMA node(s) specified by \p set
##
##  If neither ::HWLOC_MEMBIND_PROCESS nor ::HWLOC_MEMBIND_THREAD is
##  specified, the current process is assumed to be single-threaded.
##  This is the most portable form as it permits hwloc to use either
##  process-based OS functions or thread-based OS functions, depending
##  on which are available.
##
##  If ::HWLOC_MEMBIND_BYNODESET is specified, set is considered a nodeset.
##  Otherwise it's a cpuset.
##
##  \return -1 with errno set to ENOSYS if the action is not supported
##  \return -1 with errno set to EXDEV if the binding cannot be enforced
##

proc hwloc_set_membind*(topology: hwloc_topology_t; set: hwloc_const_bitmap_t;
                       policy: hwloc_membind_policy_t; flags: cint): cint {.
    importc: "hwloc_set_membind", hwloc.}
## * \brief Query the default memory binding policy and physical locality of the
##  current process or thread.
##
##  This function has two output parameters: \p set and \p policy.
##  The values returned in these parameters depend on both the \p flags
##  passed in and the current memory binding policies and nodesets in
##  the queried target.
##
##  Passing the ::HWLOC_MEMBIND_PROCESS flag specifies that the query
##  target is the current policies and nodesets for all the threads in
##  the current process.  Passing ::HWLOC_MEMBIND_THREAD specifies that
##  the query target is the current policy and nodeset for only the
##  thread invoking this function.
##
##  If neither of these flags are passed (which is the most portable
##  method), the process is assumed to be single threaded.  This allows
##  hwloc to use either process-based OS functions or thread-based OS
##  functions, depending on which are available.
##
##  ::HWLOC_MEMBIND_STRICT is only meaningful when ::HWLOC_MEMBIND_PROCESS
##  is also specified.  In this case, hwloc will check the default
##  memory policies and nodesets for all threads in the process.  If
##  they are not identical, -1 is returned and errno is set to EXDEV.
##  If they are identical, the values are returned in \p set and \p
##  policy.
##
##  Otherwise, if ::HWLOC_MEMBIND_PROCESS is specified (and
##  ::HWLOC_MEMBIND_STRICT is \em not specified), the default set
##  from each thread is logically OR'ed together.
##  If all threads' default policies are the same, \p policy is set to
##  that policy.  If they are different, \p policy is set to
##  ::HWLOC_MEMBIND_MIXED.
##
##  In the ::HWLOC_MEMBIND_THREAD case (or when neither
##  ::HWLOC_MEMBIND_PROCESS or ::HWLOC_MEMBIND_THREAD is specified), there
##  is only one set and policy; they are returned in \p set and
##  \p policy, respectively.
##
##  If ::HWLOC_MEMBIND_BYNODESET is specified, set is considered a nodeset.
##  Otherwise it's a cpuset.
##
##  If any other flags are specified, -1 is returned and errno is set
##  to EINVAL.
##

proc hwloc_get_membind*(topology: hwloc_topology_t; set: hwloc_bitmap_t;
                       policy: ptr hwloc_membind_policy_t; flags: cint): cint {.
    importc: "hwloc_get_membind", hwloc.}
## * \brief Set the default memory binding policy of the specified
##  process to prefer the NUMA node(s) specified by \p set
##
##  If ::HWLOC_MEMBIND_BYNODESET is specified, set is considered a nodeset.
##  Otherwise it's a cpuset.
##
##  \return -1 with errno set to ENOSYS if the action is not supported
##  \return -1 with errno set to EXDEV if the binding cannot be enforced
##
##  \note \p hwloc_pid_t is \p pid_t on Unix platforms,
##  and \p HANDLE on native Windows platforms.
##

proc hwloc_set_proc_membind*(topology: hwloc_topology_t; pid: hwloc_pid_t;
                            set: hwloc_const_bitmap_t;
                            policy: hwloc_membind_policy_t; flags: cint): cint {.
    importc: "hwloc_set_proc_membind", hwloc.}
## * \brief Query the default memory binding policy and physical locality of the
##  specified process.
##
##  This function has two output parameters: \p set and \p policy.
##  The values returned in these parameters depend on both the \p flags
##  passed in and the current memory binding policies and nodesets in
##  the queried target.
##
##  Passing the ::HWLOC_MEMBIND_PROCESS flag specifies that the query
##  target is the current policies and nodesets for all the threads in
##  the specified process.  If ::HWLOC_MEMBIND_PROCESS is not specified
##  (which is the most portable method), the process is assumed to be
##  single threaded.  This allows hwloc to use either process-based OS
##  functions or thread-based OS functions, depending on which are
##  available.
##
##  Note that it does not make sense to pass ::HWLOC_MEMBIND_THREAD to
##  this function.
##
##  If ::HWLOC_MEMBIND_STRICT is specified, hwloc will check the default
##  memory policies and nodesets for all threads in the specified
##  process.  If they are not identical, -1 is returned and errno is
##  set to EXDEV.  If they are identical, the values are returned in \p
##  set and \p policy.
##
##  Otherwise, \p set is set to the logical OR of all threads'
##  default set.  If all threads' default policies
##  are the same, \p policy is set to that policy.  If they are
##  different, \p policy is set to ::HWLOC_MEMBIND_MIXED.
##
##  If ::HWLOC_MEMBIND_BYNODESET is specified, set is considered a nodeset.
##  Otherwise it's a cpuset.
##
##  If any other flags are specified, -1 is returned and errno is set
##  to EINVAL.
##
##  \note \p hwloc_pid_t is \p pid_t on Unix platforms,
##  and \p HANDLE on native Windows platforms.
##

proc hwloc_get_proc_membind*(topology: hwloc_topology_t; pid: hwloc_pid_t;
                            set: hwloc_bitmap_t;
                            policy: ptr hwloc_membind_policy_t; flags: cint): cint {.
    importc: "hwloc_get_proc_membind", hwloc.}
## * \brief Bind the already-allocated memory identified by (addr, len)
##  to the NUMA node(s) specified by \p set.
##
##  If ::HWLOC_MEMBIND_BYNODESET is specified, set is considered a nodeset.
##  Otherwise it's a cpuset.
##
##  \return 0 if \p len is 0.
##  \return -1 with errno set to ENOSYS if the action is not supported
##  \return -1 with errno set to EXDEV if the binding cannot be enforced
##

proc hwloc_set_area_membind*(topology: hwloc_topology_t; `addr`: pointer; len: csize_t;
                            set: hwloc_const_bitmap_t;
                            policy: hwloc_membind_policy_t; flags: cint): cint {.
    importc: "hwloc_set_area_membind", hwloc.}
## * \brief Query the CPUs near the physical NUMA node(s) and binding policy of
##  the memory identified by (\p addr, \p len ).
##
##  This function has two output parameters: \p set and \p policy.
##  The values returned in these parameters depend on both the \p flags
##  passed in and the memory binding policies and nodesets of the pages
##  in the address range.
##
##  If ::HWLOC_MEMBIND_STRICT is specified, the target pages are first
##  checked to see if they all have the same memory binding policy and
##  nodeset.  If they do not, -1 is returned and errno is set to EXDEV.
##  If they are identical across all pages, the set and policy are
##  returned in \p set and \p policy, respectively.
##
##  If ::HWLOC_MEMBIND_STRICT is not specified, the union of all NUMA
##  node(s) containing pages in the address range is calculated.
##  If all pages in the target have the same policy, it is returned in
##  \p policy.  Otherwise, \p policy is set to ::HWLOC_MEMBIND_MIXED.
##
##  If ::HWLOC_MEMBIND_BYNODESET is specified, set is considered a nodeset.
##  Otherwise it's a cpuset.
##
##  If any other flags are specified, -1 is returned and errno is set
##  to EINVAL.
##
##  If \p len is 0, -1 is returned and errno is set to EINVAL.
##

proc hwloc_get_area_membind*(topology: hwloc_topology_t; `addr`: pointer; len: csize_t;
                            set: hwloc_bitmap_t;
                            policy: ptr hwloc_membind_policy_t; flags: cint): cint {.
    importc: "hwloc_get_area_membind", hwloc.}
## * \brief Get the NUMA nodes where memory identified by (\p addr, \p len ) is physically allocated.
##
##  Fills \p set according to the NUMA nodes where the memory area pages
##  are physically allocated. If no page is actually allocated yet,
##  \p set may be empty.
##
##  If pages spread to multiple nodes, it is not specified whether they spread
##  equitably, or whether most of them are on a single node, etc.
##
##  The operating system may move memory pages from one processor
##  to another at any time according to their binding,
##  so this function may return something that is already
##  outdated.
##
##  If ::HWLOC_MEMBIND_BYNODESET is specified in \p flags, set is
##  considered a nodeset. Otherwise it's a cpuset.
##
##  If \p len is 0, \p set is emptied.
##

proc hwloc_get_area_memlocation*(topology: hwloc_topology_t; `addr`: pointer;
                                len: csize_t; set: hwloc_bitmap_t; flags: cint): cint {.
    importc: "hwloc_get_area_memlocation", hwloc.}
## * \brief Allocate some memory
##
##  This is equivalent to malloc(), except that it tries to allocate
##  page-aligned memory from the OS.
##
##  \note The allocated memory should be freed with hwloc_free().
##

proc hwloc_alloc*(topology: hwloc_topology_t; len: csize_t): pointer {.
    importc: "hwloc_alloc", hwloc.}
## * \brief Allocate some memory on NUMA memory nodes specified by \p set
##
##  \return NULL with errno set to ENOSYS if the action is not supported
##  and ::HWLOC_MEMBIND_STRICT is given
##  \return NULL with errno set to EXDEV if the binding cannot be enforced
##  and ::HWLOC_MEMBIND_STRICT is given
##  \return NULL with errno set to ENOMEM if the memory allocation failed
##  even before trying to bind.
##
##  If ::HWLOC_MEMBIND_BYNODESET is specified, set is considered a nodeset.
##  Otherwise it's a cpuset.
##
##  \note The allocated memory should be freed with hwloc_free().
##

proc hwloc_alloc_membind*(topology: hwloc_topology_t; len: csize_t;
                         set: hwloc_const_bitmap_t;
                         policy: hwloc_membind_policy_t; flags: cint): pointer {.
    importc: "hwloc_alloc_membind", hwloc.}
## * \brief Allocate some memory on NUMA memory nodes specified by \p set
##
##  This is similar to hwloc_alloc_membind_nodeset() except that it is allowed to change
##  the current memory binding policy, thus providing more binding support, at
##  the expense of changing the current state.
##
##  If ::HWLOC_MEMBIND_BYNODESET is specified, set is considered a nodeset.
##  Otherwise it's a cpuset.
##

proc hwloc_alloc_membind_policy*(topology: hwloc_topology_t; len: csize_t;
                                set: hwloc_const_bitmap_t;
                                policy: hwloc_membind_policy_t; flags: cint): pointer {.
    importc: "hwloc_alloc_membind_policy", hwloc.}
## * \brief Free memory that was previously allocated by hwloc_alloc()
##  or hwloc_alloc_membind().
##

proc hwloc_free*(topology: hwloc_topology_t; `addr`: pointer; len: csize_t): cint {.
    importc: "hwloc_free", hwloc.}
## * @}
## * \defgroup hwlocality_setsource Changing the Source of Topology Discovery
##
##  If none of the functions below is called, the default is to detect all the objects
##  of the machine that the caller is allowed to access.
##
##  This default behavior may also be modified through environment variables
##  if the application did not modify it already.
##  Setting HWLOC_XMLFILE in the environment enforces the discovery from a XML
##  file as if hwloc_topology_set_xml() had been called.
##  Setting HWLOC_SYNTHETIC enforces a synthetic topology as if
##  hwloc_topology_set_synthetic() had been called.
##
##  Finally, HWLOC_THISSYSTEM enforces the return value of
##  hwloc_topology_is_thissystem().
##
##  @{
##
## * \brief Change which process the topology is viewed from.
##
##  On some systems, processes may have different views of the machine, for
##  instance the set of allowed CPUs. By default, hwloc exposes the view from
##  the current process. Calling hwloc_topology_set_pid() permits to make it
##  expose the topology of the machine from the point of view of another
##  process.
##
##  \note \p hwloc_pid_t is \p pid_t on Unix platforms,
##  and \p HANDLE on native Windows platforms.
##
##  \note -1 is returned and errno is set to ENOSYS on platforms that do not
##  support this feature.
##

proc hwloc_topology_set_pid*(topology: hwloc_topology_t; pid: hwloc_pid_t): cint {.
    importc: "hwloc_topology_set_pid", hwloc.}
## * \brief Enable synthetic topology.
##
##  Gather topology information from the given \p description,
##  a space-separated string of <type:number> describing
##  the object type and arity at each level.
##  All types may be omitted (space-separated string of numbers) so that
##  hwloc chooses all types according to usual topologies.
##  See also the \ref synthetic.
##
##  Setting the environment variable HWLOC_SYNTHETIC
##  may also result in this behavior.
##
##  If \p description was properly parsed and describes a valid topology
##  configuration, this function returns 0.
##  Otherwise -1 is returned and errno is set to EINVAL.
##
##  Note that this function does not actually load topology
##  information; it just tells hwloc where to load it from.  You'll
##  still need to invoke hwloc_topology_load() to actually load the
##  topology information.
##
##  \note For convenience, this backend provides empty binding hooks which just
##  return success.
##
##  \note On success, the synthetic component replaces the previously enabled
##  component (if any), but the topology is not actually modified until
##  hwloc_topology_load().
##

proc hwloc_topology_set_synthetic*(topology: hwloc_topology_t; description: cstring): cint {.
    importc: "hwloc_topology_set_synthetic", hwloc.}
## * \brief Enable XML-file based topology.
##
##  Gather topology information from the XML file given at \p xmlpath.
##  Setting the environment variable HWLOC_XMLFILE may also result in this behavior.
##  This file may have been generated earlier with hwloc_topology_export_xml() in hwloc/export.h,
##  or lstopo file.xml.
##
##  Note that this function does not actually load topology
##  information; it just tells hwloc where to load it from.  You'll
##  still need to invoke hwloc_topology_load() to actually load the
##  topology information.
##
##  \return -1 with errno set to EINVAL on failure to read the XML file.
##
##  \note See also hwloc_topology_set_userdata_import_callback()
##  for importing application-specific object userdata.
##
##  \note For convenience, this backend provides empty binding hooks which just
##  return success.  To have hwloc still actually call OS-specific hooks, the
##  ::HWLOC_TOPOLOGY_FLAG_IS_THISSYSTEM has to be set to assert that the loaded
##  file is really the underlying system.
##
##  \note On success, the XML component replaces the previously enabled
##  component (if any), but the topology is not actually modified until
##  hwloc_topology_load().
##

proc hwloc_topology_set_xml*(topology: hwloc_topology_t; xmlpath: cstring): cint {.
    importc: "hwloc_topology_set_xml", hwloc.}
## * \brief Enable XML based topology using a memory buffer (instead of
##  a file, as with hwloc_topology_set_xml()).
##
##  Gather topology information from the XML memory buffer given at \p
##  buffer and of length \p size.  This buffer may have been filled
##  earlier with hwloc_topology_export_xmlbuffer() in hwloc/export.h.
##
##  Note that this function does not actually load topology
##  information; it just tells hwloc where to load it from.  You'll
##  still need to invoke hwloc_topology_load() to actually load the
##  topology information.
##
##  \return -1 with errno set to EINVAL on failure to read the XML buffer.
##
##  \note See also hwloc_topology_set_userdata_import_callback()
##  for importing application-specific object userdata.
##
##  \note For convenience, this backend provides empty binding hooks which just
##  return success.  To have hwloc still actually call OS-specific hooks, the
##  ::HWLOC_TOPOLOGY_FLAG_IS_THISSYSTEM has to be set to assert that the loaded
##  file is really the underlying system.
##
##  \note On success, the XML component replaces the previously enabled
##  component (if any), but the topology is not actually modified until
##  hwloc_topology_load().
##

proc hwloc_topology_set_xmlbuffer*(topology: hwloc_topology_t; buffer: cstring;
                                  size: cint): cint {.
    importc: "hwloc_topology_set_xmlbuffer", hwloc.}
## * \brief Flags to be passed to hwloc_topology_set_components()
##

type
  hwloc_topology_components_flag_e* {.size: sizeof(cint).} = enum ## * \brief Blacklist the target component from being used.
                                                             ##  \hideinitializer
                                                             ##
    HWLOC_TOPOLOGY_COMPONENTS_FLAG_BLACKLIST = (1 shl 0)


## * \brief Prevent a discovery component from being used for a topology.
##
##  \p name is the name of the discovery component that should not be used
##  when loading topology \p topology. The name is a string such as "cuda".
##
##  For components with multiple phases, it may also be suffixed with the name
##  of a phase, for instance "linux:io".
##
##  \p flags should be ::HWLOC_TOPOLOGY_COMPONENTS_FLAG_BLACKLIST.
##
##  This may be used to avoid expensive parts of the discovery process.
##  For instance, CUDA-specific discovery may be expensive and unneeded
##  while generic I/O discovery could still be useful.
##

proc hwloc_topology_set_components*(topology: hwloc_topology_t; flags: culong;
                                   name: cstring): cint {.
    importc: "hwloc_topology_set_components", hwloc.}
## * @}
## * \defgroup hwlocality_configuration Topology Detection Configuration and Query
##
##  Several functions can optionally be called between hwloc_topology_init() and
##  hwloc_topology_load() to configure how the detection should be performed,
##  e.g. to ignore some objects types, define a synthetic topology, etc.
##
##  @{
##
## * \brief Flags to be set onto a topology context before load.
##
##  Flags should be given to hwloc_topology_set_flags().
##  They may also be returned by hwloc_topology_get_flags().
##

type
  hwloc_topology_flags_e* {.size: sizeof(cint).} = enum ## * \brief Detect the whole system, ignore reservations, include disallowed objects.
                                                   ##
                                                   ##  Gather all resources, even if some were disabled by the administrator.
                                                   ##  For instance, ignore Linux Cgroup/Cpusets and gather all processors and memory nodes.
                                                   ##
                                                   ##  When this flag is not set, PUs and NUMA nodes that are disallowed are not added to the topology.
                                                   ##  Parent objects (package, core, cache, etc.) are added only if some of their children are allowed.
                                                   ##  All existing PUs and NUMA nodes in the topology are allowed.
                                                   ##
                                                   ## hwloc_topology_get_allowed_cpuset() and
                                                   ## hwloc_topology_get_allowed_nodeset()
                                                   ##  are equal to the root object cpuset and nodeset.
                                                   ##
                                                   ##  When this flag is set, the actual sets of allowed PUs and NUMA nodes are given
                                                   ##  by
                                                   ## hwloc_topology_get_allowed_cpuset() and
                                                   ## hwloc_topology_get_allowed_nodeset().
                                                   ##  They may be smaller than the root object cpuset and nodeset.
                                                   ##
                                                   ##  If the current topology is exported to XML and reimported later, this flag
                                                   ##  should be set again in the reimported topology so that disallowed resources
                                                   ##  are reimported as well.
                                                   ##  \hideinitializer
                                                   ##
    HWLOC_TOPOLOGY_FLAG_INCLUDE_DISALLOWED = (1 shl 0), ## * \brief Assume that the selected backend provides the topology for the
                                                   ##  system on which we are running.
                                                   ##
                                                   ##  This forces
                                                   ## hwloc_topology_is_thissystem() to return 1, i.e. makes hwloc assume that
                                                   ##  the selected backend provides the topology for the system on which we are running,
                                                   ##  even if it is not the OS-specific backend but the XML backend for instance.
                                                   ##  This means making the binding functions actually call the OS-specific
                                                   ##  system calls and really do binding, while the XML backend would otherwise
                                                   ##  provide empty hooks just returning success.
                                                   ##
                                                   ##  Setting the environment variable HWLOC_THISSYSTEM may also result in the
                                                   ##  same behavior.
                                                   ##
                                                   ##  This can be used for efficiency reasons to first detect the topology once,
                                                   ##  save it to an XML file, and quickly reload it later through the XML
                                                   ##  backend, but still having binding functions actually do bind.
                                                   ##  \hideinitializer
                                                   ##
    HWLOC_TOPOLOGY_FLAG_IS_THISSYSTEM = (1 shl 1), ## * \brief Get the set of allowed resources from the local operating system even if the topology was loaded from XML or synthetic description.
                                              ##
                                              ##  If the topology was loaded from XML or from a synthetic string,
                                              ##  restrict it by applying the current process restrictions such as
                                              ##  Linux Cgroup/Cpuset.
                                              ##
                                              ##  This is useful when the topology is not loaded directly from
                                              ##  the local machine (e.g. for performance reason) and it comes
                                              ##  with all resources, while the running process is restricted
                                              ##  to only parts of the machine.
                                              ##
                                              ##  This flag is ignored unless
                                              ## ::HWLOC_TOPOLOGY_FLAG_IS_THISSYSTEM is
                                              ##  also set since the loaded topology must match the underlying machine
                                              ##  where restrictions will be gathered from.
                                              ##
                                              ##  Setting the environment variable HWLOC_THISSYSTEM_ALLOWED_RESOURCES
                                              ##  would result in the same behavior.
                                              ##  \hideinitializer
                                              ##
    HWLOC_TOPOLOGY_FLAG_THISSYSTEM_ALLOWED_RESOURCES = (1 shl 2)


## * \brief Set OR'ed flags to non-yet-loaded topology.
##
##  Set a OR'ed set of ::hwloc_topology_flags_e onto a topology that was not yet loaded.
##
##  If this function is called multiple times, the last invokation will erase
##  and replace the set of flags that was previously set.
##
##  The flags set in a topology may be retrieved with hwloc_topology_get_flags()
##

proc hwloc_topology_set_flags*(topology: hwloc_topology_t; flags: culong): cint {.
    importc: "hwloc_topology_set_flags", hwloc.}
## * \brief Get OR'ed flags of a topology.
##
##  Get the OR'ed set of ::hwloc_topology_flags_e of a topology.
##
##  \return the flags previously set with hwloc_topology_set_flags().
##

proc hwloc_topology_get_flags*(topology: hwloc_topology_t): culong {.
    importc: "hwloc_topology_get_flags", hwloc.}
## * \brief Does the topology context come from this system?
##
##  \return 1 if this topology context was built using the system
##  running this program.
##  \return 0 instead (for instance if using another file-system root,
##  a XML topology file, or a synthetic topology).
##

proc hwloc_topology_is_thissystem*(topology: hwloc_topology_t): cint {.
    importc: "hwloc_topology_is_thissystem", hwloc.}
## * \brief Flags describing actual discovery support for this topology.

type
  hwloc_topology_discovery_support* {.bycopy.} = object
    pu*: uint8                ## * \brief Detecting the number of PU objects is supported.
    ## * \brief Detecting the number of NUMA nodes is supported.
    numa*: uint8              ## * \brief Detecting the amount of memory in NUMA nodes is supported.
    numa_memory*: uint8       ## * \brief Detecting and identifying PU objects that are not available to the current process is supported.
    disallowed_pu*: uint8     ## * \brief Detecting and identifying NUMA nodes that are not available to the current process is supported.
    disallowed_numa*: uint8


## * \brief Flags describing actual PU binding support for this topology.
##
##  A flag may be set even if the feature isn't supported in all cases
##  (e.g. binding to random sets of non-contiguous objects).
##

type
  hwloc_topology_cpubind_support* {.bycopy.} = object
    set_thisproc_cpubind*: uint8 ## * Binding the whole current process is supported.
    ## * Getting the binding of the whole current process is supported.
    get_thisproc_cpubind*: uint8 ## * Binding a whole given process is supported.
    set_proc_cpubind*: uint8  ## * Getting the binding of a whole given process is supported.
    get_proc_cpubind*: uint8  ## * Binding the current thread only is supported.
    set_thisthread_cpubind*: uint8 ## * Getting the binding of the current thread only is supported.
    get_thisthread_cpubind*: uint8 ## * Binding a given thread only is supported.
    set_thread_cpubind*: uint8 ## * Getting the binding of a given thread only is supported.
    get_thread_cpubind*: uint8 ## * Getting the last processors where the whole current process ran is supported
    get_thisproc_last_cpu_location*: uint8 ## * Getting the last processors where a whole process ran is supported
    get_proc_last_cpu_location*: uint8 ## * Getting the last processors where the current thread ran is supported
    get_thisthread_last_cpu_location*: uint8


## * \brief Flags describing actual memory binding support for this topology.
##
##  A flag may be set even if the feature isn't supported in all cases
##  (e.g. binding to random sets of non-contiguous objects).
##

type
  hwloc_topology_membind_support* {.bycopy.} = object
    set_thisproc_membind*: uint8 ## * Binding the whole current process is supported.
    ## * Getting the binding of the whole current process is supported.
    get_thisproc_membind*: uint8 ## * Binding a whole given process is supported.
    set_proc_membind*: uint8  ## * Getting the binding of a whole given process is supported.
    get_proc_membind*: uint8  ## * Binding the current thread only is supported.
    set_thisthread_membind*: uint8 ## * Getting the binding of the current thread only is supported.
    get_thisthread_membind*: uint8 ## * Binding a given memory area is supported.
    set_area_membind*: uint8  ## * Getting the binding of a given memory area is supported.
    get_area_membind*: uint8  ## * Allocating a bound memory area is supported.
    alloc_membind*: uint8     ## * First-touch policy is supported.
    firsttouch_membind*: uint8 ## * Bind policy is supported.
    bind_membind*: uint8      ## * Interleave policy is supported.
    interleave_membind*: uint8 ## * Next-touch migration policy is supported.
    nexttouch_membind*: uint8 ## * Migration flags is supported.
    migrate_membind*: uint8   ## * Getting the last NUMA nodes where a memory area was allocated is supported
    get_area_memlocation*: uint8


## * \brief Set of flags describing actual support for this topology.
##
##  This is retrieved with hwloc_topology_get_support() and will be valid until
##  the topology object is destroyed.  Note: the values are correct only after
##  discovery.
##

type
  hwloc_topology_support* {.bycopy.} = object
    discovery*: ptr hwloc_topology_discovery_support
    cpubind*: ptr hwloc_topology_cpubind_support
    membind*: ptr hwloc_topology_membind_support


## * \brief Retrieve the topology support.
##
##  Each flag indicates whether a feature is supported.
##  If set to 0, the feature is not supported.
##  If set to 1, the feature is supported, but the corresponding
##  call may still fail in some corner cases.
##
##  These features are also listed by hwloc-info \--support
##

proc hwloc_topology_get_support*(topology: hwloc_topology_t): ptr hwloc_topology_support {.
    importc: "hwloc_topology_get_support", hwloc.}
## * \brief Type filtering flags.
##
##  By default, most objects are kept (::HWLOC_TYPE_FILTER_KEEP_ALL).
##  Instruction caches, I/O and Misc objects are ignored by default (::HWLOC_TYPE_FILTER_KEEP_NONE).
##  Die and Group levels are ignored unless they bring structure (::HWLOC_TYPE_FILTER_KEEP_STRUCTURE).
##
##  Note that group objects are also ignored individually (without the entire level)
##  when they do not bring structure.
##

type
  hwloc_type_filter_e* {.size: sizeof(cint).} = enum ## * \brief Keep all objects of this type.
                                                ##
                                                ##  Cannot be set for ::HWLOC_OBJ_GROUP (groups are designed only to add more structure to the topology).
                                                ##  \hideinitializer
                                                ##
    HWLOC_TYPE_FILTER_KEEP_ALL = 0, ## * \brief Ignore all objects of this type.
                                 ##
                                 ##  The bottom-level type ::HWLOC_OBJ_PU, the ::HWLOC_OBJ_NUMANODE type, and
                                 ##  the top-level type ::HWLOC_OBJ_MACHINE may not be ignored.
                                 ##  \hideinitializer
                                 ##
    HWLOC_TYPE_FILTER_KEEP_NONE = 1, ## * \brief Only ignore objects if their entire level does not bring any structure.
                                  ##
                                  ##  Keep the entire level of objects if at least one of these objects adds
                                  ##  structure to the topology. An object brings structure when it has multiple
                                  ##  children and it is not the only child of its parent.
                                  ##
                                  ##  If all objects in the level are the only child of their parent, and if none
                                  ##  of them has multiple children, the entire level is removed.
                                  ##
                                  ##  Cannot be set for I/O and Misc objects since the topology structure does not matter there.
                                  ##  \hideinitializer
                                  ##
    HWLOC_TYPE_FILTER_KEEP_STRUCTURE = 2, ## * \brief Only keep likely-important objects of the given type.
                                       ##
                                       ##  It is only useful for I/O object types.
                                       ##  For ::HWLOC_OBJ_PCI_DEVICE and ::HWLOC_OBJ_OS_DEVICE, it means that only objects
                                       ##  of major/common kinds are kept (storage, network, OpenFabrics, Intel MICs, CUDA,
                                       ##  OpenCL, NVML, and displays).
                                       ##  Also, only OS devices directly attached on PCI (e.g. no USB) are reported.
                                       ##  For ::HWLOC_OBJ_BRIDGE, it means that bridges are kept only if they have children.
                                       ##
                                       ##  This flag equivalent to ::HWLOC_TYPE_FILTER_KEEP_ALL for Normal, Memory and Misc types
                                       ##  since they are likely important.
                                       ##  \hideinitializer
                                       ##
    HWLOC_TYPE_FILTER_KEEP_IMPORTANT = 3


## * \brief Set the filtering for the given object type.
##

proc hwloc_topology_set_type_filter*(topology: hwloc_topology_t;
                                    `type`: hwloc_obj_type_t;
                                    filter: hwloc_type_filter_e): cint {.
    importc: "hwloc_topology_set_type_filter", hwloc.}
## * \brief Get the current filtering for the given object type.
##

proc hwloc_topology_get_type_filter*(topology: hwloc_topology_t;
                                    `type`: hwloc_obj_type_t;
                                    filter: ptr hwloc_type_filter_e): cint {.
    importc: "hwloc_topology_get_type_filter", hwloc.}
## * \brief Set the filtering for all object types.
##
##  If some types do not support this filtering, they are silently ignored.
##

proc hwloc_topology_set_all_types_filter*(topology: hwloc_topology_t;
    filter: hwloc_type_filter_e): cint {.importc: "hwloc_topology_set_all_types_filter",
                                      hwloc.}
## * \brief Set the filtering for all CPU cache object types.
##
##  Memory-side caches are not involved since they are not CPU caches.
##

proc hwloc_topology_set_cache_types_filter*(topology: hwloc_topology_t;
    filter: hwloc_type_filter_e): cint {.importc: "hwloc_topology_set_cache_types_filter",
                                      hwloc.}
## * \brief Set the filtering for all CPU instruction cache object types.
##
##  Memory-side caches are not involved since they are not CPU caches.
##

proc hwloc_topology_set_icache_types_filter*(topology: hwloc_topology_t;
    filter: hwloc_type_filter_e): cint {.importc: "hwloc_topology_set_icache_types_filter",
                                      hwloc.}
## * \brief Set the filtering for all I/O object types.
##

proc hwloc_topology_set_io_types_filter*(topology: hwloc_topology_t;
                                        filter: hwloc_type_filter_e): cint {.
    importc: "hwloc_topology_set_io_types_filter", hwloc.}
## * \brief Set the topology-specific userdata pointer.
##
##  Each topology may store one application-given private data pointer.
##  It is initialized to \c NULL.
##  hwloc will never modify it.
##
##  Use it as you wish, after hwloc_topology_init() and until hwloc_topolog_destroy().
##
##  This pointer is not exported to XML.
##

proc hwloc_topology_set_userdata*(topology: hwloc_topology_t; userdata: pointer) {.
    importc: "hwloc_topology_set_userdata", hwloc.}
## * \brief Retrieve the topology-specific userdata pointer.
##
##  Retrieve the application-given private data pointer that was
##  previously set with hwloc_topology_set_userdata().
##

proc hwloc_topology_get_userdata*(topology: hwloc_topology_t): pointer {.
    importc: "hwloc_topology_get_userdata", hwloc.}
## * @}
## * \defgroup hwlocality_tinker Modifying a loaded Topology
##  @{
##
## * \brief Flags to be given to hwloc_topology_restrict().

type
  hwloc_restrict_flags_e* {.size: sizeof(cint).} = enum ## * \brief Remove all objects that became CPU-less.
                                                   ##  By default, only objects that contain no PU and no memory are removed.
                                                   ##  \hideinitializer
                                                   ##
    HWLOC_RESTRICT_FLAG_REMOVE_CPULESS = (1 shl 0), ## * \brief Restrict by nodeset instead of CPU set.
                                               ##  Only keep objects whose nodeset is included or partially included in the given set.
                                               ##  This flag may not be used with ::HWLOC_RESTRICT_FLAG_BYNODESET.
                                               ##
    HWLOC_RESTRICT_FLAG_ADAPT_MISC = (1 shl 1), ## * \brief Move I/O objects to ancestors if their parents are removed during restriction.
                                           ##  If this flag is not set, I/O devices and bridges are removed when their parents are removed.
                                           ##  \hideinitializer
                                           ##
    HWLOC_RESTRICT_FLAG_ADAPT_IO = (1 shl 2),
    HWLOC_RESTRICT_FLAG_BYNODESET = (1 shl 3), ## * \brief Remove all objects that became Memory-less.
                                          ##  By default, only objects that contain no PU and no memory are removed.
                                          ##  This flag may only be used with ::HWLOC_RESTRICT_FLAG_BYNODESET.
                                          ##  \hideinitializer
                                          ##
    HWLOC_RESTRICT_FLAG_REMOVE_MEMLESS = (1 shl 4), ## * \brief Move Misc objects to ancestors if their parents are removed during restriction.
                                               ##  If this flag is not set, Misc objects are removed when their parents are removed.
                                               ##  \hideinitializer
                                               ##


## * \brief Restrict the topology to the given CPU set or nodeset.
##
##  Topology \p topology is modified so as to remove all objects that
##  are not included (or partially included) in the CPU set \p set.
##  All objects CPU and node sets are restricted accordingly.
##
##  If ::HWLOC_RESTRICT_FLAG_BYNODESET is passed in \p flags,
##  \p set is considered a nodeset instead of a CPU set.
##
##  \p flags is a OR'ed set of ::hwloc_restrict_flags_e.
##
##  \note This call may not be reverted by restricting back to a larger
##  set. Once dropped during restriction, objects may not be brought
##  back, except by loading another topology with hwloc_topology_load().
##
##  \return 0 on success.
##
##  \return -1 with errno set to EINVAL if the input set is invalid.
##  The topology is not modified in this case.
##
##  \return -1 with errno set to ENOMEM on failure to allocate internal data.
##  The topology is reinitialized in this case. It should be either
##  destroyed with hwloc_topology_destroy() or configured and loaded again.
##

proc hwloc_topology_restrict*(topology: hwloc_topology_t;
                             set: hwloc_const_bitmap_t; flags: culong): cint {.
    importc: "hwloc_topology_restrict", hwloc.}
## * \brief Flags to be given to hwloc_topology_allow().

type
  hwloc_allow_flags_e* {.size: sizeof(cint).} = enum ## * \brief Mark all objects as allowed in the topology.
                                                ##
                                                ##  \p cpuset and \p nođeset given to hwloc_topology_allow() must be \c NULL.
                                                ##  \hideinitializer
    HWLOC_ALLOW_FLAG_ALL = (1 shl 0), ## * \brief Only allow objects that are available to the current process.
                                 ##
                                 ##  The topology must have ::HWLOC_TOPOLOGY_FLAG_IS_THISSYSTEM so that the set
                                 ##  of available resources can actually be retrieved from the operating system.
                                 ##
                                 ##  \p cpuset and \p nođeset given to hwloc_topology_allow() must be \c NULL.
                                 ##  \hideinitializer
    HWLOC_ALLOW_FLAG_LOCAL_RESTRICTIONS = (1 shl 1), ## * \brief Allow a custom set of objects, given to hwloc_topology_allow() as \p cpuset and/or \p nodeset parameters.
                                                ##  \hideinitializer
    HWLOC_ALLOW_FLAG_CUSTOM = (1 shl 2)


## * \brief Change the sets of allowed PUs and NUMA nodes in the topology.
##
##  This function only works if the ::HWLOC_TOPOLOGY_FLAG_INCLUDE_DISALLOWED
##  was set on the topology. It does not modify any object, it only changes
##  the sets returned by hwloc_topology_get_allowed_cpuset() and
##  hwloc_topology_get_allowed_nodeset().
##
##  It is notably useful when importing a topology from another process
##  running in a different Linux Cgroup.
##
##  \p flags must be set to one flag among ::hwloc_allow_flags_e.
##
##  \note Removing objects from a topology should rather be performed with
##  hwloc_topology_restrict().
##

proc hwloc_topology_allow*(topology: hwloc_topology_t;
                          cpuset: hwloc_const_cpuset_t;
                          nodeset: hwloc_const_nodeset_t; flags: culong): cint {.
    importc: "hwloc_topology_allow", hwloc.}
## * \brief Add a MISC object as a leaf of the topology
##
##  A new MISC object will be created and inserted into the topology at the
##  position given by parent. It is appended to the list of existing Misc children,
##  without ever adding any intermediate hierarchy level. This is useful for
##  annotating the topology without actually changing the hierarchy.
##
##  \p name is supposed to be unique across all Misc objects in the topology.
##  It will be duplicated to setup the new object attributes.
##
##  The new leaf object will not have any \p cpuset.
##
##  \return the newly-created object
##
##  \return \c NULL on error.
##
##  \return \c NULL if Misc objects are filtered-out of the topology (::HWLOC_TYPE_FILTER_KEEP_NONE).
##
##  \note If \p name contains some non-printable characters, they will
##  be dropped when exporting to XML, see hwloc_topology_export_xml() in hwloc/export.h.
##

proc hwloc_topology_insert_misc_object*(topology: hwloc_topology_t;
                                       parent: hwloc_obj_t; name: cstring): hwloc_obj_t {.
    importc: "hwloc_topology_insert_misc_object", hwloc.}
## * \brief Allocate a Group object to insert later with hwloc_topology_insert_group_object().
##
##  This function returns a new Group object.
##  The caller should (at least) initialize its sets before inserting the object.
##  See hwloc_topology_insert_group_object().
##
##  The \p subtype object attribute may be set to display something else
##  than "Group" as the type name for this object in lstopo.
##  Custom name/value info pairs may be added with hwloc_obj_add_info() after
##  insertion.
##
##  The \p kind group attribute should be 0. The \p subkind group attribute may
##  be set to identify multiple Groups of the same level.
##
##  It is recommended not to set any other object attribute before insertion,
##  since the Group may get discarded during insertion.
##
##  The object will be destroyed if passed to hwloc_topology_insert_group_object()
##  without any set defined.
##

proc hwloc_topology_alloc_group_object*(topology: hwloc_topology_t): hwloc_obj_t {.
    importc: "hwloc_topology_alloc_group_object", hwloc.}
## * \brief Add more structure to the topology by adding an intermediate Group
##
##  The caller should first allocate a new Group object with hwloc_topology_alloc_group_object().
##  Then it must setup at least one of its CPU or node sets to specify
##  the final location of the Group in the topology.
##  Then the object can be passed to this function for actual insertion in the topology.
##
##  The group \p dont_merge attribute may be set to prevent the core from
##  ever merging this object with another object hierarchically-identical.
##
##  Either the cpuset or nodeset field (or both, if compatible) must be set
##  to a non-empty bitmap. The complete_cpuset or complete_nodeset may be set
##  instead if inserting with respect to the complete topology
##  (including disallowed, offline or unknown objects).
##
##  It grouping several objects, hwloc_obj_add_other_obj_sets() is an easy way
##  to build the Group sets iteratively.
##
##  These sets cannot be larger than the current topology, or they would get
##  restricted silently.
##
##  The core will setup the other sets after actual insertion.
##
##  \return The inserted object if it was properly inserted.
##
##  \return An existing object if the Group was discarded because the topology already
##  contained an object at the same location (the Group did not add any locality information).
##  Any name/info key pair set before inserting is appended to the existing object.
##
##  \return \c NULL if the insertion failed because of conflicting sets in topology tree.
##
##  \return \c NULL if Group objects are filtered-out of the topology (::HWLOC_TYPE_FILTER_KEEP_NONE).
##
##  \return \c NULL if the object was discarded because no set was initialized in the Group
##  before insert, or all of them were empty.
##

proc hwloc_topology_insert_group_object*(topology: hwloc_topology_t;
                                        group: hwloc_obj_t): hwloc_obj_t {.
    importc: "hwloc_topology_insert_group_object", hwloc.}
## * \brief Setup object cpusets/nodesets by OR'ing another object's sets.
##
##  For each defined cpuset or nodeset in \p src, allocate the corresponding set
##  in \p dst and add \p src to it by OR'ing sets.
##
##  This function is convenient between hwloc_topology_alloc_group_object()
##  and hwloc_topology_insert_group_object(). It builds the sets of the new Group
##  that will be inserted as a new intermediate parent of several objects.
##

proc hwloc_obj_add_other_obj_sets*(dst: hwloc_obj_t; src: hwloc_obj_t): cint {.
    importc: "hwloc_obj_add_other_obj_sets", hwloc.}
## * @}

##  high-level helpers

##  inline code of some functions above

##  exporting to XML or synthetic

##  distances

##  topology diffs

##  deprecated headers
## * \defgroup hwlocality_helper_topology_sets CPU and node sets of entire topologies
##  @{
##
## * \brief Get complete CPU set
##
##  \return the complete CPU set of logical processors of the system.
##
##  \note The returned cpuset is not newly allocated and should thus not be
##  changed or freed; hwloc_bitmap_dup() must be used to obtain a local copy.
##
##  \note This is equivalent to retrieving the root object complete CPU-set.
##

proc hwloc_topology_get_complete_cpuset*(topology: hwloc_topology_t): hwloc_const_cpuset_t {.
    importc: "hwloc_topology_get_complete_cpuset", hwloc.}
## * \brief Get topology CPU set
##
##  \return the CPU set of logical processors of the system for which hwloc
##  provides topology information. This is equivalent to the cpuset of the
##  system object.
##
##  \note The returned cpuset is not newly allocated and should thus not be
##  changed or freed; hwloc_bitmap_dup() must be used to obtain a local copy.
##
##  \note This is equivalent to retrieving the root object CPU-set.
##

proc hwloc_topology_get_topology_cpuset*(topology: hwloc_topology_t): hwloc_const_cpuset_t {.
    importc: "hwloc_topology_get_topology_cpuset", hwloc.}
## * \brief Get allowed CPU set
##
##  \return the CPU set of allowed logical processors of the system.
##
##  \note If the topology flag ::HWLOC_TOPOLOGY_FLAG_INCLUDE_DISALLOWED was not set,
##  this is identical to hwloc_topology_get_topology_cpuset(), which means
##  all PUs are allowed.
##
##  \note If ::HWLOC_TOPOLOGY_FLAG_INCLUDE_DISALLOWED was set, applying
##  hwloc_bitmap_intersects() on the result of this function and on an object
##  cpuset checks whether there are allowed PUs inside that object.
##  Applying hwloc_bitmap_and() returns the list of these allowed PUs.
##
##  \note The returned cpuset is not newly allocated and should thus not be
##  changed or freed, hwloc_bitmap_dup() must be used to obtain a local copy.
##

proc hwloc_topology_get_allowed_cpuset*(topology: hwloc_topology_t): hwloc_const_cpuset_t {.
    importc: "hwloc_topology_get_allowed_cpuset", hwloc.}
## * \brief Get complete node set
##
##  \return the complete node set of memory of the system.
##
##  \note The returned nodeset is not newly allocated and should thus not be
##  changed or freed; hwloc_bitmap_dup() must be used to obtain a local copy.
##
##  \note This is equivalent to retrieving the root object complete nodeset.
##

proc hwloc_topology_get_complete_nodeset*(topology: hwloc_topology_t): hwloc_const_nodeset_t {.
    importc: "hwloc_topology_get_complete_nodeset", hwloc.}
## * \brief Get topology node set
##
##  \return the node set of memory of the system for which hwloc
##  provides topology information. This is equivalent to the nodeset of the
##  system object.
##
##  \note The returned nodeset is not newly allocated and should thus not be
##  changed or freed; hwloc_bitmap_dup() must be used to obtain a local copy.
##
##  \note This is equivalent to retrieving the root object nodeset.
##

proc hwloc_topology_get_topology_nodeset*(topology: hwloc_topology_t): hwloc_const_nodeset_t {.
    importc: "hwloc_topology_get_topology_nodeset", hwloc.}
## * \brief Get allowed node set
##
##  \return the node set of allowed memory of the system.
##
##  \note If the topology flag ::HWLOC_TOPOLOGY_FLAG_INCLUDE_DISALLOWED was not set,
##  this is identical to hwloc_topology_get_topology_nodeset(), which means
##  all NUMA nodes are allowed.
##
##  \note If ::HWLOC_TOPOLOGY_FLAG_INCLUDE_DISALLOWED was set, applying
##  hwloc_bitmap_intersects() on the result of this function and on an object
##  nodeset checks whether there are allowed NUMA nodes inside that object.
##  Applying hwloc_bitmap_and() returns the list of these allowed NUMA nodes.
##
##  \note The returned nodeset is not newly allocated and should thus not be
##  changed or freed, hwloc_bitmap_dup() must be used to obtain a local copy.
##

proc hwloc_topology_get_allowed_nodeset*(topology: hwloc_topology_t): hwloc_const_nodeset_t {.
    importc: "hwloc_topology_get_allowed_nodeset", hwloc.}
## * @}
