#ifndef GLIDE_MASK_INC

#define GLIDE_MASK_OCEAN           -2
#define GLIDE_MASK_LAND            -1

!Set if ice is below the ice dynamics thickness limit
#define GLIDE_MASK_THIN_ICE        1

!Grounded Ice
#define GLIDE_MASK_INTERIOR        2

!*sfp* true for lateral boundary of any kind
#define GLIDE_MASK_BOUNDARY        4

#define GLIDE_MASK_SHELF           8
#define GLIDE_MASK_STREAM          16
#define GLIDE_MASK_GROUNDING_LINE  32
#define GLIDE_MASK_STREAM_MARGIN   64
#define GLIDE_MASK_LAND_MARGIN     128  
#define GLIDE_MASK_MARINE_EDGE     256

!Derived mask values
#define GLIDE_MASK_SHELF_FRONT ior(GLIDE_MASK_MARINE_EDGE, GLIDE_MASK_SHELF)

#define GLIDE_IS_OCEAN(mask)       (mask .eq. GLIDE_MASK_OCEAN)
#define GLIDE_IS_LAND(mask)        (mask .eq. GLIDE_MASK_LAND)
#define GLIDE_HAS_ICE(mask)        (mask .gt. 0)
#define GLIDE_IS_THIN(mask)        (iand(mask,GLIDE_MASK_THIN_ICE) .gt. 0 .and. GLIDE_HAS_ICE(mask))
#define GLIDE_IS_FLOAT(mask)       (iand(mask,GLIDE_MASK_SHELF) .gt. 0 .and. GLIDE_HAS_ICE(mask))

#define GLIDE_IS_GROUND(mask)      (iand(mask,GLIDE_MASK_INTERIOR) .gt. 0 .and. GLIDE_HAS_ICE(mask))

!True if at an ice boundary that is next to ocean and if the ice is floating
#define GLIDE_IS_CALVING(mask)     (iand(mask,GLIDE_MASK_SHELF_FRONT) .eq. GLIDE_MASK_SHELF_FRONT .and. GLIDE_HAS_ICE(mask))

!True if at an ice boundary that is next to ocean, regardless of whether the ice is floating or not
#define GLIDE_IS_MARINE_ICE_EDGE(mask) (iand(mask,GLIDE_MASK_MARINE_EDGE) .gt. 0 .and. GLIDE_HAS_ICE(mask))

#define GLIDE_IS_GROUNDING_LINE(mask) (iand(mask, GLIDE_MASK_GROUNDING_LINE) .gt. 0 .and. GLIDE_HAS_ICE(mask))

#endif

