#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  COMPONENTS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        TABULATED
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK3
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        1
#define  USER_DEF_PARAMETERS            8

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO
#define  PRINT_TO_FILE                  YES
#define  INTERNAL_BOUNDARY              YES

/* -- user-defined parameters (labels) -- */

#define  CHI                            0
#define  ALPHA                          1
#define  TCL                            2
#define  MACH                           3
#define  NCL                            4
#define  X0                             5
#define  Y0                             6
#define  Z0                             7

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY                   CONST_mp
#define  UNIT_LENGTH                    CONST_pc
#define  UNIT_VELOCITY                  1.e5

/* [End] user-defined constants (do not change this line) */
