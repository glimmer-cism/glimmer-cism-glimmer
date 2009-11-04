module plume_global

  implicit none

  !parameter used to declare double-precision reals
  integer,parameter::kdp = kind(1.0d0)
  integer,parameter::ksp = kind(1.0e0)

  ! set array sizes to maximum (data storage actually used set later)
  ! lamb = maximum number of vertical levels for ambient fluid t,s,rho 
  ! lice = maximum no of frazil size classes
  ! linf = maximum no of inflows
  ! ldr = maximum number of regions in which drag coefficient is varied
  ! lsm = maximum number of regions in which smoothing is varied

  integer lm,ln
  integer lamb, lice, linf,ldr,lsm,ldec
  integer lxin,lyin

  parameter (lm=1000,ln=1000)
  parameter (lamb=302,lice=10,linf=9,ldr=3,lsm=3,ldec=26)
  parameter (lxin=1000,lyin=1000)

  !variables for cell dimensions

  real(kind=kdp) hx,hy

  real(kind=kdp) dx(lm),dxu(lm)
  real(kind=kdp) rdx(lm),rdxu(lm)   ! reciprocols of dx and dxu

  real(kind=kdp) dy(ln),dyv(ln)
  real(kind=kdp) rdy(ln),rdyv(ln)   ! reciprocols of dy and dyv
  ! variables for plume thickness, interface position, 
  ! ice shelf bottom depth and position

  !plume thickness
  real(kind=kdp) pdep(lm,ln)

  !interface position
  real(kind=kdp) ipos(lm,ln)

  !ice shelf bottom depth
  real(kind=kdp) bpos(lm,ln)

  !cell incices demarking status of cell

  !     jc-land/sea
  !     jcd-wet/dry
  !      _fl      =>  newly-wet
  !      _negdep  =>  negative depths
  !      _fseed   =>  frazil seeded


  integer jc(lm,ln),jcd(lm,ln),jcd_u(lm,ln),jcd_v(lm,ln)
  integer jcd_fl(lm,ln),jcd_negdep(lm,ln),jcd_fseed(lm,ln)

  !ice-related variables

  real(kind=kdp) c(lm,ln,lice),ca(lm,ln,lice),ctot(lm,ln)
  real(kind=kdp) ctota(lm,ln),tf(lm,ln),bmelt(lm,ln)
  real(kind=kdp) btemp(lm,ln),bsalt(lm,ln),ctempd(lm,ln)

  !ice-related parameters

  real(kind=kdp) lat,c0,ci,nu0,pr,sc,fta,ftb,ftc,ti,si,nus,kt,ks
  real(kind=kdp) ar,eps,nbar,r(lice),re(lice),thi(lice),vol(lice)
  real(kind=kdp) wi(lice),cmin(lice),cseed(lice),nuss(lice)

  !frazil interaction terms

  real(kind=kdp) fmelt(lm,ln,lice),fppn(lm,ln,lice)
  real(kind=kdp) fnuc(lm,ln,lice)

  !inflow-related variables and initial thickness

  logical inflag(linf)
  integer infloa,infloe,knfloa,knfloe,intrin(lm,ln)

  real(kind=kdp) saltinf(lm,ln),tempinf(lm,ln),depinf(lm,ln)
  real(kind=kdp) intr(lm,ln,linf),intra(lm,ln,linf)
  real(kind=kdp) meltinf,cinf(lice),cinftot,depinffix,depinit

  !array limits actually used
  
  integer m,n,namb,nice,ninfmin,ninfmax

  !switches
  
  logical mixlayer,restart,nonlin,horturb,entrain,basmelt,frazil
  logical rholinear,thermobar,intrace,vardrag,topedit,tangle,negfrz
  integer entype

  !general parameters

  real(kind=kdp) pi,dcr,g,dt,gdt,fdt,small,edepth,mdepth,fdepth,septol
  real(kind=kdp) ah,kh,dzincr,temptop,tempbot,saltbot,salttop
  real(kind=kdp) tgrad,sgrad,wcdep,gldep,ifdep,rho0,rhoi
  real(kind=kdp) cdb,cdbvar,ef,cl

  !scalar fields

  real(kind=kdp) rhop(lm,ln),temp(lm,ln),tempa(lm,ln),tins(lm,ln)
  real(kind=kdp) salt(lm,ln),salta(lm,ln),rhoamb(lm,ln)
  real(kind=kdp) samb(lamb),tamb(lamb),rhovf(lamb),entr(lm,ln)
  real(kind=kdp) atemp(lm,ln),asalt(lm,ln),drag(lm,ln)
  logical drflag(ldr)

  !topography parameters

  character(len=6) context
  logical smflag(lsm)
  integer bathtype,kcorn,rad,bsmoothit,smoothit(lsm)
  real(kind=kdp) cweight,nweight

  !eddy viscosities

  real(kind=kdp) ahdx(lm),ahdxu(lm),ahdy(ln),ahdyv(ln)

  !transport and velocities

  real(kind=kdp) u(lm,ln),ua(lm,ln),v(lm,ln),va(lm,ln)
  real(kind=kdp) su(lm,ln),sv(lm,ln),u0(lm,ln),v0(lm,ln)
  real(kind=kdp) u0a(lm,ln),v0a(lm,ln),tang(lm,ln)
  
end module plume_global

