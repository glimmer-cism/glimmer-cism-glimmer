#ifndef GLIDE_MASK_INC

#define GLIDE_MASK_OCEAN           -2
#define GLIDE_MASK_LAND            -1
#define GLIDE_MASK_BOUNDARY        0

!Set if ice is below the ice dynamics thickness limit
#define GLIDE_MASK_THIN_ICE        1

!Grounded Ice
#define GLIDE_MASK_INTERIOR        2

#define GLIDE_MASK_SHELF           4
#define GLIDE_MASK_STREAM          8
#define GLIDE_MASK_GROUNDING_LINE  16
#define GLIDE_MASK_STREAM_MARGIN   32
#define GLIDE_MASK_LAND_MARGIN     64
#define GLIDE_MASK_SHELF_FRONT     128
#define GLIDE_MASK_MARINE_EDGE     256


#define GLIDE_IS_GROUNDING_LINE(mask) (iand(mask, GLIDE_MASK_GROUNDING_LINE) .gt. 0 .and. GLIDE_HAS_ICE(mask))
#define GLIDE_IS_OCEAN(mask)       (mask .eq. GLIDE_MASK_OCEAN)
#define GLIDE_IS_LAND(mask)        (mask .eq. GLIDE_MASK_LAND)
#define GLIDE_HAS_ICE(mask)        (mask .gt. 0)
#define GLIDE_IS_THIN(mask)        (iand(mask,GLIDE_MASK_THIN_ICE) .gt. 0 .and. GLIDE_HAS_ICE(mask))
#define GLIDE_IS_FLOAT(mask)       (iand(mask,GLIDE_MASK_SHELF) .gt. 0 .and. GLIDE_HAS_ICE(mask))
#define GLIDE_IS_GROUND(mask)      (iand(mask,GLIDE_MASK_INTERIOR) .gt. 0 .and. GLIDE_HAS_ICE(mask))
#define GLIDE_IS_CALVING(mask)     (iand(mask,GLIDE_MASK_SHELF_FRONT) .gt. 0 .and. GLIDE_HAS_ICE(mask))
#define GLIDE_IS_MARINE_ICE_EDGE(mask) (iand(mask,GLIDE_MASK_MARINE_EDGE) .gt. 0 .and. GLIDE_HAS_ICE(mask))
#define GLIDE_IS_SHELF_FRONT(mask) (iand(mask,GLIDE_MASK_SHELF_FRONT) .gt. 0 .and. GLIDE_HAS_ICE(mask))
#endif

