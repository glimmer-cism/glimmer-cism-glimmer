#ifndef GLIDE_MASK_INC

!Bits 1:0 - Ice presence (0 if no ice)
!NOTE: If bit 2 is activated, bit 1 must be activated
#define GLIDE_MASK_HAS_ICE             1
#define GLIDE_MASK_THIN_ICE            3
#define GLIDE_ICE_PRESENCE_BITS        3

!Bits 4:2 - Type of base (Land or ocean - grounding line has both bits on)
!The 16 bit specifies grounding line - I am setting this up so that those points are treated as grounded
!as well.
#define GLIDE_MASK_LAND                4
#define GLIDE_MASK_OCEAN               8
#define GLIDE_MASK_GROUNDING_LINE      16 
#define GLIDE_BASE_TYPE_BITS           28

!Bits 6:5 - Presence of a boundary condition.  If bit 4 is activated, it's a domain boundary or margin.
!The meaning of bit 5 depends on whether bit 4 is activated.  If it is not, bit 5 determines whether
!a dirichelt boundary was specified.  If it is, bit 5 is on if the boundary is on the domain interior
!(i.e. if it is a land margin or calving front)
!If this is a margin, the type of margin is dependent on whether the ice is on land or ocean.
!A value of 0 for these bits means that we are on the interior.
#define GLIDE_MASK_DOMAIN_BOUNDARY     32
#define GLIDE_MASK_MARGIN              96
#define GLIDE_MASK_DIRICHLET_BND       64
#define GLIDE_BOUNDARY_TYPE_BITS       96

!Checks for an iceless square

!Checks for open ocean with no ice
#define GLIDE_IS_OCEAN(mask)       (iand(mask, GLIDE_MASK_OCEAN) == GLIDE_MASK_OCEAN .and. GLIDE_NO_ICE(mask))

!Checks for land with no ice
#define GLIDE_IS_LAND(mask)        (iand(mask, GLIDE_MASK_LAND) == GLIDE_MASK_LAND .and. GLIDE_NO_ICE(mask))

!Checks for the presence of any ice, dynamic or not
#define GLIDE_HAS_ICE(mask)        (iand(mask, GLIDE_MASK_HAS_ICE) == GLIDE_MASK_HAS_ICE)

!Checks for a lack of ice
#define GLIDE_NO_ICE(mask)         (iand(mask, GLIDE_MASK_HAS_ICE) == 0)

!Checks for the presence of ice that is below the ice dynamics limit
#define GLIDE_IS_THIN(mask)          (iand(mask,GLIDE_MASK_THIN_ICE)  == GLIDE_MASK_THIN_ICE)

!Checks for any ice, dynamic or not, that is on an ice shelf.
#define GLIDE_IS_FLOAT(mask)         (iand(mask,GLIDE_MASK_OCEAN) == GLIDE_MASK_OCEAN .and. GLIDE_HAS_ICE(mask))

!Checks for any ice, dynamic or not, that is grounded
#define GLIDE_IS_GROUND(mask)      (iand(mask,GLIDE_MASK_LAND) == GLIDE_MASK_LAND .and. GLIDE_HAS_ICE(mask))

!Checks for any ice, dynamic or not, that is on the grounding line
#define GLIDE_IS_GROUNDING_LINE(mask) (iand(mask, GLIDE_MASK_GROUNDING_LINE) == GLIDE_MASK_GROUNDING_LINE .and. GLIDE_HAS_ICE(mask))

!Checks for any ice, dynamic or not, that is either floating *or* on the grounding line
#define GLIDE_IS_FLOAT_OR_GNDLINE(mask) (iand(mask, 24) > 0 .and. GLIDE_HAS_ICE(mask))

!Checks for any ice, dynamic or not, that is either grounded *or* on the grounding line
#define GLIDE_IS_GROUND_OR_GNDLINE(mask) (iand(mask, 20) > 0 .and. GLIDE_HAS_ICE(mask))

!Checks whether this is a domain boundary that requires upwinding (either a margin or a computational boundary)
#define GLIDE_IS_UPWINDED_BOUNDARY(mask) (iand(mask, GLIDE_MASK_DOMAIN_BOUNDARY) == GLIDE_MASK_DOMAIN_BOUNDARY)

!Checks whether this is an ice margin (thickness jumps from 0 to non-zero at this point)
#define GLIDE_IS_MARGIN(mask) (iand(mask, GLIDE_MASK_MARGIN) == GLIDE_MASK_MARGIN)

!Checks whether this is a margin in contact with the ocean, floating or not
#define GLIDE_IS_MARINE_ICE_EDGE(mask)    (GLIDE_IS_MARGIN(mask) .and. GLIDE_IS_FLOAT_OR_GNDLINE(mask))

!Checks whether this is a margin in contact with the ocean
#define GLIDE_IS_CALVING(mask)            (GLIDE_IS_MARGIN(mask) .and. GLIDE_IS_FLOAT(mask))

!Checks whether this is a land margin
#define GLIDE_IS_LAND_MARGIN(mask)        (GLIDE_IS_MARGIN(mask) .and. GLIDE_IS_GOUNDED(mask))


#endif

