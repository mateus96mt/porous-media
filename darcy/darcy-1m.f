c      program stabilzed mixed and hybrid fem
c      Lagrange        
c
c     ************************************************************
c     *                                                          *
c     *              * * *     DARCY FLOW      * * *             *
c     *                                                          *
c     *                                                          *
c     *       DUAL HIBRID STABILIZED FINITE ELEMENT METHODS      *
c     *                                                          *
c     *                 LOCALLY CONSERVATIVE                     *
c     *                                                          *
c     *                                                          *
c     *      IURY IGREJA, CRISTIANE FARIA E ABIMAEL LOULA        *
c     *                                                          *
c     *                     April 2013                           *
c     ************************************************************
c
c.... program to set storage capacity, precision and input/output units
c
c
c     PROBLEMA DE DARCY COM 2 MULTIPLICADORES 
c
c     (Possibilidade de testes usando interpolacoes de diferentes graus para u e p)
c
c     - 2*µ Æu + grad p = f em K
c                eps p + div u  = 0
c                             u = g em dK
c
c
c      Escolha dos Multiplicadores 
c                  _
c      \lambda_u = u 
c
c      betah = beta0/h 
c
c.... program to set storage capacity, precision and input/output units
c
      common /bpoint/ mfirst,mlast,ilast,mtot,iprec
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi

      character*4 ia 
      parameter (ndim=400000000) 
      common a(ndim) 
      common /dictn/ ia(10000000)
c
c        mfirst = address of first available word in blank common
c        mlast  = address of last available word in blank common
c        mtot   = total storage allocated to blank common
c        iprec  = precision flag; eq.1, single precision
c                                 eq.2, double precision
c
      iin= 8
      iout= 11
      iecho=12
      ioupp=13
      iout1=14
c      itest2=15
c      ierrb0=16
c      ierrb=17
c      ierrbi=18

c               LEITURA
      open(unit=iin, file= 'darcy-1m - 40e.dat',status='old')
      open(unit=iecho, file= 'darcy.eco')	
      open(unit=iout, file= 'errofem-le.con')      
      open(unit=ioupp, file= 'erro.local.proj.con')
      open(unit=iout1, file= 'erro.interpolante.con')      
c      open(unit=itest2, file= 'teste2.con')
c arquivos de estudo de convergencia: log(erro)x beta0
c      open(unit=ierrb0, file= 'errobetalp.con')      
c      open(unit=ierrb, file= 'errobeta.con')      
c      open(unit=ierrbi, file= 'errobetai.con')      

c
      mfirst = 1 
      ilast  = 0 
      mlast  = ndim 
      mtot   = ndim 
      iprec  = 2 
c
c        iin    = input unit number
c        iecho  = output unit of input data
c        ioupp  = output unit of post-processed displacements
c
c
       call lpgm
c
c.... system-dependent unit/file specifications
c
c
         close(iin)
         close(iecho)
	 close(iout)
	 close(ioupp)
	 close(iout1)
c	 close(itest2)
c         close(ierrb0)
c         close(ierrb)
c         close(ierrbi)

c
         write(*,*)'fim'
      stop
      end
c
**********************************************************************
      subroutine lpgm
c
c.... LPGM - a linear static finite element analysis program for 
c            Petrov Galerkin methods : global driver
c
      real*8 zero,pt1667,pt25,pt5,one,two,three,four,five,six,
     & tempf
      character*4 title
c
c.... remove above card for single-precision operation
c
c.... catalog of common statements
c
      common /bpoint/ mfirst,mlast,ilast,mtot,iprec
      common /colhtc/ neq
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /etimec/ etime(6)
      common /genelc/ n,nel(3),incel(3),inc(3)                          
      common /genflc/ tempf(6,20),nf,numgpf,nincf(3),incf(3)
      common /info  / iexec,iprtin,irank,nsd,numnp,ndof,nlvect,
     &                numeg,nmultp,nedge
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
      common /labels/ labeld(3),label1(16),label2(3)
      common /spoint/ mpd,mpx,mpid,mpf,mpdiag,mpngrp,
     &                mpalhs,mpbrhs,mped,index,
     &                mpdlhs

      common /titlec/ title(20)
      character*4 ia 
      common a(1) 
      common /dictn/ ia(1) 

c
c.... input phase
c
      call echo
c
  100 continue
      do 200 i=1,6
  200 etime(i) = 0.0
      read(iin,1000) title
      if (title(1).eq.'*end') return
      read(iin,2000) iexec,iprtin,irank,
     &           nsd,numnp,ndof,nlvect,numeg,nedge,npar
c
      write(iecho,3000) title , iexec,iprtin
      write(iecho,4000) irank , nsd, numnp,  ndof,
     &                  nlvect,numeg,nedge,npar
c
       nmultp = nedge*npar
c
c.... initialization phase
c
c....    set memory pointers for static data arrays,
c        and call associated input routines 
c
      mpd    = mpoint('d       ',ndof  ,nmultp ,0     ,iprec)
      mpx    = mpoint('x       ',nsd   ,numnp  ,0     ,iprec)
      mped   = mpoint('ideg    ',2*ndof  ,nedge  ,0     ,1)
      mpid   = mpoint('id      ',ndof  ,nmultp ,0     ,1)
c
      if (nlvect.eq.0) then
         mpf = 1
      else
         mpf = mpoint('f       ',ndof  ,nmultp ,nlvect,iprec)
      endif
c
c.... input coordinate data
c
      call coord(a(mpx),nsd,numnp,iprtin)
c
c.... input boundary condition data and establish equation numbers
c
      call bcedge(a(mped),a(mpid),npar,nedge,ndof,nmultp,neq,iprtin)
c
c.... input nodal force and prescribed kinematic boundary-value data
c
      if (nlvect.gt.0) call input(a(mpf),ndof,nmultp,0,nlvect,
     &                            iprtin)
c
c.... allocate memory for idiag array and clear 
c
      mpdiag = mpoint('idiag   ',neq   ,0     ,0     ,1)
      call iclear(a(mpdiag),neq)
c
      mpngrp = mpoint('ngrp    ',numeg ,0     ,0     ,1)
c
c.... input element data
c
      call elemnt('input___',a(mpngrp))
c
c.... determine addresses of diagonals in left-hand-side matrix
c
      call diag(a(mpdiag),neq,nalhs)
c
c.... allocate memory for global equation system
c
      mpalhs = mpoint('alhs    ',nalhs,0,0,iprec)
      mpdlhs = mpoint('dlhs    ',nalhs,0,0,iprec)
      mpbrhs = mpoint('brhs    ',neq  ,0,0,iprec)
      meanbw = nalhs/neq
      nwords = mtot - mlast + mfirst - 1
c
c.... write equation system data
c
      write(iecho,*)nwords,mtot,mlast,mfirst
      write(iecho,5000) title,neq,nalhs,meanbw,nwords
c
c.... solution phase
c
      if (iexec.eq.1) call driver(neq,nalhs)
c
c.... print memory-pointer dictionary
c
c      call prtdc
c
c
c.... print elapsed time summary
c
      go to 100
c
 1000 format(20a4)
 2000 format(17i10)
 3000 format(///,20a4///
     &' e x e c u t i o n   c o n t r o l   i n f o r m a t i o n '//5x,
     &' execution code  . . . . . . . . . . . . . . (iexec ) = ',i10//5x,
     &'    eq. 0, data check                                   ',   /5x,
     &'    eq. 1, execution                                    ',  //5x,
     &' input data print code . . . . . . . . . . . (iprtin) = ',i10//5x,
     &'    eq. 0, print nodal and element input data           ',   /5x,
     &'    eq. 1, do not print nodal and element input data    ',   /5x)
 4000 format(5x,
     &' rank check code . . . . . . . . . . . . . . (irank ) = ',i10//5x,
     &'    eq. 0, do not perform rank check                    ',   /5x,
     &'    eq. 1, print numbers of zero and nonpositive pivots ',   /5x,
     &'    eq. 2, print all pivots                             ',  //5x,
     &' number of space dimensions  . . . . . . . .(nsd   ) = ',i10//5x,
     &' number of nodal points  . . . . . . . . . .(numnp ) = ',i10//5x,
     &' number of nodal degrees-of-freedom  . . . .(ndof  ) = ',i10//5x,
     &' number of load vectors  . . . . . . . . . .(nlvect) = ',i10//5x,
     &' number of element groups  . . . . . . . . .(numeg ) = ',i10//5x,
     &' number of edge. . . .  .  . . . . . . . . .(nedge ) = ',i10//5x,
     &' number of nodal point of multiply . . . . .(npar  ) = ',i10//5x)
4500  format(///,20a4//5x,
     &'nsd = ',i10,5x,'nmultp = ',i10,5x,'ndof = ',i10,5x,
     &'nlvect = ',i10,5x,'numeg = ',i10)
 5000 format(///,20a4///
     &' e q u a t i o n    s y s t e m    d a t a              ',  //5x,
     &' number of equations . . . . . . . . . . . . (neq   ) = ',i8//5x,   
     &' number of terms in left-hand-side matrix  . (nalhs ) = ',i8//5x,
     &' mean half bandwidth . . . . . . . . . . . . (meanbw) = ',i8//5x,
     &' total length of blank common required . . . (nwords) = ',i8    )
c
      end
c***********************************************************************
c**** new **********************************************************************
      subroutine driver(neq,nalhs)
c
c.... solution driver program 
c
      common /etimec/ etime(6)
      common /info  / iexec,iprtin,irank,nsd,numnp,ndof,nlvect,
     &                numeg,nmultp,nedge

      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb  
     & ierrbi
      common /spoint/ mpd,mpx,mpid,mpf,mpdiag,mpngrp,
     &                mpalhs,mpbrhs,mped,index,
     &                mpdlhs
      character*4 ia 
      common a(1) 
      common /dictn/ ia(1)
c
c      clear left and right hand side
c
      call clear(a(mpalhs),nalhs)
      call clear(a(mpdlhs),nalhs)
c
      call clear(a(mpbrhs),neq)
c
c
      call elemnt('form_stb',a(mpngrp))
c
c      account the nodal forces in the r.h.s.
c
      if (nlvect.gt.0)
     &   call load(a(mpid),a(mpf),a(mpbrhs),ndof,nmultp,nlvect)
c
c	clear displacement array
c
      call clear(a(mpd),ndof*nmultp)
c
      if (nlvect.gt.0)
     &   call ftod(a(mpid),a(mpd),a(mpf),ndof,nmultp,nlvect)
c 
c
c      form the l.h.s and r.h.s. at element level
c
      call elemnt('form_lrs',a(mpngrp))
c
c      factorization of the stiffnes matrix
c
      if(neq.eq.0) go to 1111
      call factor(a(mpalhs),a(mpdiag),neq)
c      call factns(a(mpalhs),a(mpdlhs),a(mpdiag),neq)
c
c      back substitution
c
      call back(a(mpalhs),a(mpbrhs),a(mpdiag),neq)
c      call backns(a(mpalhs),a(mpdlhs),a(mpbrhs),a(mpdiag),neq)
c
 1111 continue  
c
      call btod(a(mpid),a(mpd),a(mpbrhs),ndof,nmultp)
c
c.... write output 
c
c
         call printd(' M u l t i p l y                          ',
     &               a(mpd),ndof,nmultp,iecho)
c         write(166,*)'Multiply'
c         call   printdX(' M u l t i p l y                          ',
c     &               a(mpd),ndof,nmultp,166)
c
c
c
c	post-processing phase
c
      call elemnt('pos_proc',a(mpngrp))
c
c
c
  100 continue
      return
      end 

c***********************************************************************
c             scopo do programa
c***********************************************************************

      subroutine elemnt(task,ngrp)
c
c.... program to calculate element task number
c
      character*8 task,eltask(4) 
      dimension ngrp(*)
      common /info  / iexec,iprtin,irank,nsd,numnp,ndof,nlvect,
     &                numeg,nmultp,nedge

      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
      character*4 ia 
      common na(1) 
      common /dictn/ ia(1) 
      data ntask,    eltask
     &    /    4,'input___',
     &           'form_stb',
     &           'form_lrs',
     &           'pos_proc'/
c
      do 100 i=1,ntask
      if (task.eq.eltask(i)) itask = i
  100 continue
c
      do 200 neg=1,numeg
c
      if (itask.eq.1) then
         mpnpar = mpoint('npar    ',16   ,0,0,1)
         ngrp(neg) = mpnpar
         call elcard(na(mpnpar),neg)
      else
         mpnpar = ngrp(neg)
      endif
c
      ntype  = na(mpnpar)
      call elmlib(ntype,mpnpar,itask,neg)
  200 continue
c
      return
      end
c**** new **********************************************************************
      subroutine elcard(npar,neg)
c
c.... program to read element group control card
c
      dimension npar(*)
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
      read(iin,1000) (npar(i),i=1,16)
      write(iecho,2000) neg
c
      return
c
 1000 format(16i10)
 2000 format(//,' e l e m e n t   g r o u p   d a t a         ',  //5x,
     &' element group number  . . . . . . . . . . (neg   ) = ',i10/// )
c
      end
c**** new **********************************************************************
      subroutine elmlib(ntype,mpnpar,itask,neg)
c
c.... program to call element routines
c
      common a(1)
c
      go to (10,20) ntype
c
c    POISSON PROBLEM - KINEMATIC FORMULATION
c
  10  continue
      call darcy(itask,a(mpnpar),a(mpnpar+16),neg)
  20  continue
      return
      end
c**** new **********************************************************************
      subroutine darcy(itask,npar,mp,neg)
c___________________________________________________________
c
c..... program to set storage and call tasks for the 
c             primal mixed Poisson  problem
c     with continuous temperature and discontinuous flux
c___________________________________________________________
      dimension npar(*),mp(*)
      common /bpoint/ mfirst,mlast,ilast,mtot,iprec
      common /colhtc/ neq      
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
      common /info  / iexec,iprtin,irank,nsd,numnp,ndof,nlvect,
     &                numeg,nmultp,nedge

      common /spoint/ mpd,mpx,mpid,mpf,mpdiag,mpngrp,
     &                mpalhs,mpbrhs,mped,index,
     &                mpdlhs
      character*4 ia
      common a(1)
      common /dictn/ ia(1)
c
c
      mw     = 1
      mdet   = 2
      mshl   = 3
      mshg   = 4
c
      mc     = 5
      mgrav  = 6
      mien   = 7
      mmat   = 8
      mlm    = 9
c
      mxl    = 10
      mdl    = 11
c
c	pointers for condensation
c
      mipar  = 12
      mlado  = 13
c
      melefd = 14
      melred = 15
      mdlf   = 17	  
      mdlp   = 18
c
c    inetrais nas arestas
c
      mdetpn = 19
      mshlpn = 20
      mshgpn = 21
      mwpn   = 22
c
c.....pointers for boundary itegrals
c
      mdside = 23
      mxls   = 24
      midlsd = 25
c
c    geoemetria das arstas
c      
      mdetn  = 26
      mshln  = 27
      mshgn  = 28
      mwn    = 29
c
c    valores na fronteira
c
      mdetb  = 30
      mshlb  = 31
      mshgb  = 32
      mddis  = 33
      mdsfl  = 34
c
      mdetp  = 35
      mshlp  = 36
      mshgp  = 37
      mwp    = 38

c
      mdetc  = 39
      mshlc  = 40
      mshgc  = 41
      mwc    = 42
	  
c
c    matrizes e vetores da formulação hibrida
c
      melmbb  = 43
      melmcb  = 44
      melmcc  = 45
c
      melfbb  = 46
      melfcb  = 47
      melfcc  = 48
	  
c
c
      mshlpsd = 49
      mshgpsd = 50
      mshlcsd = 51
      mshgcsd = 52
c
      melmdb  = 53
c termo não simétrico      
      melmbc  = 54
c
      mshlsd = 55
      mshgsd = 56
c
cc      mepres = 57
c
      mmatside = 57
      mmcont = 58
c
      ntype  = npar( 1)
      numel  = npar( 2) !numero de elementos globais da malha   
      numat  = npar( 3) !numero de grupo de materiais 
      nint   = npar( 4) !numeros de pontos de integracao do elemento 
      nen    = npar( 5) !Numero de nós do elemento geométrico
      nencon = npar( 6) !Numero de nós do elemento para o parametro u
      nenp   = npar( 7) !Numero de nós do elemento para o parametro p
      npars  = npar( 8) !Numero de nós do multiplicador em cada lado do elemento (nosso caso o w)
      nints  = npar( 9) !numero de pontos de integracao em cada aresta (unidimensional)
      nface  = npar(10) !numero de aresta do elemento (como o nside) só que esta variável é usada
c                       !para gerar integração reduzida no termo penalizado.  
c

      write(iecho,*)'ntype  =',ntype
      write(iecho,*)'numel  =',numel
      write(iecho,*)'numat  =',numat
      write(iecho,*)'nint   =',nint
      write(iecho,*)'nen    =',nen
      write(iecho,*)'nencon =',nencon  
      write(iecho,*)'nenp   =',nenp
      write(iecho,*)'npars  =',npars
      write(iecho,*)'nints  =',nints
      write(iecho,*)'nface  =',nface
c nen = número de nós da geometria (nen=3 triangle, nen=4 quadrado)
      if((nen.eq.3).or.(nen.eq.4))then
        nenlad = 2 !número de nós do lado do elemento (geometrico)
        if(nen.eq.3) nside = 3
        if(nen.eq.4) nside = 4
      else
        stop
      endif
c elementos quadrangulares
c nnods = numero de pontos no lado para a pressao
      if(nenp.eq.1)  nnods = 1  
      if(nenp.eq.4)  nnods = 2  
      if(nenp.eq.8)  nnods = 3
      if(nenp.eq.9)  nnods = 3      
      if(nenp.eq.16) nnods = 4
      if(nenp.eq.25) nnods = 5
      if(nenp.eq.36) nnods = 6
      if(nenp.eq.49) nnods = 7
      if(nenp.eq.64) nnods = 8
c
c    elementos triangulares
c
      if(nenp.eq.3)  nnods  = 2
      if(nenp.eq.6)  nnods  = 3
      if(nenp.eq.10) nnods  = 4
      if(nenp.eq.15) nnods  = 5
      if(nenp.eq.21) nnods  = 6
      if(nenp.eq.28) nnods  = 7      
c
c nnodc = no de ptos no lado para a velocidade
c
      if(nencon.eq.1)  nnodc = 1
      if(nencon.eq.4)  nnodc = 2
      if(nencon.eq.8)  nnodc = 3
      if(nencon.eq.9)  nnodc = 3
      if(nencon.eq.16) nnodc = 4
      if(nencon.eq.25) nnodc = 5
      if(nencon.eq.36) nnodc = 6
      if(nencon.eq.49) nnodc = 7
      if(nencon.eq.64) nnodc = 8
c
      if(nencon.eq.3)  nnodc = 2
      if(nencon.eq.6)  nnodc = 3
      if(nencon.eq.10) nnodc = 4
      if(nencon.eq.15) nnodc = 5
      if(nencon.eq.21) nnodc = 6
      if(nencon.eq.28) nnodc = 7
c
c
      nintb=nside*nints
c
c      parameters for post processing
c
	  nodsp=npars*nside
c
c
c.... set element parameters
c
      ndimc  = 9
      ned    = 1  ! Numero de graus de liberdade pressao
      ncon   = 2  ! Numero de graus de liberdade veolicade
      nee    = nodsp*ndof  ! Está associado ao multiplicador
      neep   = nenp*ned    ! está associado a p 
      necon  = nencon*ncon ! esta associado a u
	  nepcon = neep+necon  ! associado a (u,p)
      neesq  = nee*nee 
      nesd   = 2
      nrowsh = 3
      ngrav  = 6
c
       if (itask.eq.1) then
c
c....... set memory pointers
c
c
c        note:  the mp array is stored directly after the npar array,
c               beginning at location mpnpar + 16 of blank common.
c               the variable "junk" is not used subsequently.
c
         junk       = mpoint('mp      ',80     ,0     ,0     ,1)
c
         mp(mw    ) = mpoint('w       ',nint   ,0     ,0     ,iprec)
         mp(mdet  ) = mpoint('det     ',nint   ,0     ,0     ,iprec)
         mp(mshl  ) = mpoint('shl     ',nrowsh ,nen   ,nint  ,iprec)
         mp(mshg  ) = mpoint('shg     ',nrowsh ,nen   ,nint  ,iprec)
         mp(mc    ) = mpoint('c       ',ndimc  ,numat ,0     ,iprec)
         mp(mgrav ) = mpoint('grav    ',ngrav  ,0     ,0     ,iprec)
c
         mp(mien  ) = mpoint('ien     ',nen    ,numel ,0     ,1)
         mp(mmat  ) = mpoint('mat     ',numel  ,0     ,0     ,1)
c
         mp(mlm   ) = mpoint('lm      ',ndof   ,nodsp ,numel ,1)
         mp(mxl   ) = mpoint('xl      ',nesd   ,nen   ,0     ,iprec)
         mp(mdl   ) = mpoint('dl      ',ndof   ,nodsp ,0     ,iprec)
c
         mp(mipar ) = mpoint('ipar    ',nodsp  ,numel ,0,     1)
         mp(mlado)  = mpoint('lado    ',nside  ,numel ,0     ,1)
c
c	
         mp(melefd) = mpoint('eleffd  ',nee    ,nee   ,0     ,iprec)
         mp(melred) = mpoint('elresd  ',nee    ,0     ,0     ,iprec)
         mp(mdlf  ) = mpoint('dlf     ',ncon   ,nencon,0     ,iprec)
         mp(mdlp  ) = mpoint('dlp     ',ned    ,nenp  ,0     ,iprec)
c
         mp(mdetpn) = mpoint('detpn   ',nints  ,0     ,0     ,iprec)
         mp(mshlpn) = mpoint('shlpn   ',nrowsh ,npars ,nints ,iprec)
         mp(mshgpn) = mpoint('shgpn   ',nrowsh ,npars ,nints ,iprec)
         mp(mwpn  ) = mpoint('wpn     ',nints  ,0     ,0     ,iprec)
c
         mp(mdside) = mpoint('idside  ',nside  ,nnods ,0     ,1    )
         mp(mxls  ) = mpoint('xls     ',nesd   ,nenlad,0     ,iprec)
         mp(midlsd) = mpoint('idlsd   ',nnods  ,0     ,0     ,1    )
         mp(mdetn ) = mpoint('detn    ',nints  ,0     ,0     ,iprec)
         mp(mshln ) = mpoint('shln    ',3      ,nnods ,nints ,iprec)
         mp(mshgn ) = mpoint('shgn    ',3      ,nnods ,nints ,iprec)
         mp(mwn   ) = mpoint('wn      ',nints  ,0     ,0     ,iprec)
c
         mp(mdetb ) = mpoint('detb    ',nints  ,0     ,0     ,iprec)
         mp(mshlb ) = mpoint('shlb    ',3      ,nenlad,nints ,iprec)
         mp(mshgb ) = mpoint('shgb    ',3      ,nenlad,nints ,iprec)
c
         mp(mddis ) = mpoint('ddis    ',ned    ,nenp  ,numel ,iprec)
         mp(mdsfl ) = mpoint('dsfl    ',ncon   ,nencon,numel ,iprec)
c
         mp(mwp   ) = mpoint('wp      ',nint   ,0     ,0     ,iprec)
         mp(mdetp ) = mpoint('detp    ',nint   ,0     ,0     ,iprec)
         mp(mshlp ) = mpoint('shlp    ',nrowsh ,nenp  ,nint  ,iprec)
         mp(mshgp ) = mpoint('shgp    ',nrowsh ,nenp  ,nint  ,iprec)
c
         mp(mwc   ) = mpoint('wc      ',nint   ,0     ,0     ,iprec)
         mp(mdetc ) = mpoint('detc    ',nint   ,0     ,0     ,iprec)
         mp(mshlc ) = mpoint('shlc    ',nrowsh ,nencon,nint  ,iprec)
         mp(mshgc ) = mpoint('shgc    ',nrowsh ,nencon,nint  ,iprec)
c
c
         mp(melmbb) = mpoint('elmbb   ',nepcon ,nepcon,0     ,iprec)
         mp(melmcb) = mpoint('elmcb   ',nepcon ,nee   ,0     ,iprec)
         mp(melmcc) = mpoint('elmcc   ',nee    ,nee   ,0     ,iprec)
c
         mp(melfbb) = mpoint('elfbb   ',nepcon ,0     ,0     ,iprec)
         mp(melfcb) = mpoint('elfcb   ',nee    ,0     ,0     ,iprec)
         mp(melfcc) = mpoint('elfcc   ',nee    ,0     ,0     ,iprec)
c
c
         mp(mshlpsd) = mpoint('shlpsd  ',nrowsh ,nenp  ,nintb ,iprec)
         mp(mshgpsd) = mpoint('shgpsd  ',nrowsh ,nenp  ,nintb ,iprec)
         mp(mshlcsd) = mpoint('shlcsd  ',nrowsh ,nencon,nintb ,iprec)
         mp(mshgcsd) = mpoint('shgcsd  ',nrowsh ,nencon,nintb ,iprec)
c
         mp(melmdb) = mpoint('elmdb   ',nee    ,nepcon ,0     ,iprec)
         mp(melmbc) = mpoint('elmbc   ',nee    ,nepcon ,0     ,iprec)
c
         mp(mshlsd) = mpoint('shlsd   ',nrowsh ,nen   ,nintb ,iprec)
         mp(mshgsd) = mpoint('shgsd   ',nrowsh ,nen   ,nintb ,iprec)
c
cc         mp(mepres) =   mpoint('epres   ',nee+2 ,numel   ,0    ,iprec)
c
         mp(mmatside) = mpoint('matside ',nedge  ,0     ,0     ,1)
         mp(mmcont)   = mpoint('mcont   ',nedge  ,0     ,0     ,1)
c
      endif
c
c.... task calls
c
      if (itask.gt.4) return
      go to (100,200,300,400),itask
c
  100 continue
c
c.... input element data ('input___')
c
      call darcy1(a(mp(mshl  )),a(mp(mw  )),
     &           a(mp(mc    )),a(mp(mgrav )),
     &           a(mp(mien  )),a(mp(mmat  )),
     &           a(mp(mmatside)),a(mp(mmcont)),
     &           a(mpid      ),a(mp(mlm   )),
     &           a(mpdiag    ),a(mp(mipar )),
     &           a(mpx       ),a(mp(mlado)),
     &           a(mp(mshlc )),a(mp(mwc   )),
     &           a(mp(mshlpn)),a(mp(mwpn  )),
     &           a(mp(mshln )),a(mp(mwn   )),
     &           a(mp(mshlb )),a(mp(mshlp )),
     &           a(mp(mwp   )),a(mp(mdside)),
c     
     &           a(mp(mshlpsd)),a(mp(mshlcsd)),
     &           a(mp(mshlsd)),
c     
     &           ntype ,numel ,numat ,
     &           nint  ,nrowsh,nesd  ,
     &           nen   ,ndof  ,ned   ,
     &           iprtin,numnp ,ncon  ,
     &           nencon,necon ,nints ,
     &           nnods ,nnodc ,nenlad,npars ,
     &           nenp  ,nside ,nodsp ,nedge )
c
      return
c
  200 continue
c
c     ('form_stb')
c
      read(iin,*) index,iwrite
c
c               comentar
      write(iecho,*)'index = ',index,'iwrite = ',iwrite
c      call darcy0(a(mp(mien  )),a(mpx       ),a(mp(mxl )),
c     &           a(mpd       ),a(mp(mdl   )),a(mp(mmat  )),
c     &           a(mp(mmatside)),a(mp(mmcont)),
c     &           a(mp(mdet  )),a(mp(mshl  )),a(mp(mshg  )),
c     &           a(mp(mw    )),a(mp(mc    )),
c     &           a(mp(mgrav )),a(mp(mipar )),a(mp(mlado )),
c     &           a(mp(mdetc )),a(mp(mshlc )),a(mp(mshgc )),
c     &           a(mp(melefd)),a(mp(melred)),a(mp(mshln )),
c     &           a(mp(mshgn )),a(mp(mwn   )),a(mp(mdetn )),
c     &           a(mp(mdetb )),a(mp(mshlb )),a(mp(mshgb )),
c     &           a(mp(mdetpn)),a(mp(mshlpn)),a(mp(mshgpn)),
c     &           a(mp(mdside)),a(mp(mxls  )),a(mp(midlsd)),
c     &           a(mp(mdsfl )),a(mp(mddis )),a(mp(mdetp )),
c     &           a(mp(mshlp )),a(mp(mshgp )),
     
c     &           a(mp(melmbb)),
c     &           a(mp(melmcb)),
c     &           a(mp(melfbb)),a(mp(melfcb)),
c     &           a(mp(melmdb)),
c
c     &           a(mped      ),
c
c     &           a(mp(mshlpsd)),a(mp(mshlcsd)),
c     &           a(mp(mshgpsd)),a(mp(mshgcsd)),
c          
c     &           a(mp(mshlsd)),
c     &           a(mp(mshgsd)),
c
c     &           numel ,neesq ,nen   ,
c     &           nsd   ,nesd  ,nint  ,
c     &           neg   ,nrowsh,ned   ,
c     &           nee   ,numnp ,ndof  ,
c     &           ncon  ,nencon,necon ,
c     &           neep  ,nints ,nnods ,
c     &           nenlad,npars ,nside  ,
c     &           nenp  ,nodsp ,index ,
c     &   	     nepcon,nface )
c

c               comentar

c      call flnorm(a(mp(mien )),a(mpx       ),a(mp(mxl   )),
c     &            a(mpd      ),a(mp(mdl   )),a(mp(mmat  )),
c     &            a(mp(mmatside)),
c     &            a(mp(mc   )),a(mp(mipar )),a(mp(mdlf )) ,
c     &            a(mp(mdlp )),a(mp(mdsfl )),a(mp(mdet  )),
c     &            a(mp(mshl )),a(mp(mshg  )),a(mp(mw    )),
c     &            a(mp(mdetc )),a(mp(mshlc)),a(mp(mshgc )),
c     &            a(mp(mddis )),a(mp(mdetp)),a(mp(mshlp )),
c     &            a(mp(mshgp )),a(mp(mlado)),
c
c     &            a(mped      ),
c
c     &            a(mp(mshln )),a(mp(mshgn )),
c     &            a(mp(mdetn )),a(mp(mshlb )),a(mp(mshgb )),
c     &            a(mp(mdetpn)),a(mp(mshlpn)),a(mp(mshgpn)),
c     &            a(mp(mdside)),a(mp(mxls  )),a(mp(midlsd)),
c     &            a(mp(mgrav )),a(mp(mwn   )),
c
c     &            numel ,neesq ,nen   ,nsd   ,
c     &            nesd  ,nint  ,neg   ,nrowsh,
c     &            ned   ,nee   ,numnp ,ndof  ,
c     &            ncon  ,nencon,necon ,index ,
c     &            nints ,iwrite,ioupp ,ierrb0,
c     &            nenp  ,nside ,nnods ,nenlad,
c     &            npars ,nmultp,nodsp)
c
c      return
c
 300  continue
c
c.... form element stiffnes matrix and force vector and
c     assemble them into global left and right hand-side
c     ('form_lrs')
c
      call darcy2(a(mp(mien  )),a(mpx       ),a(mp(mxl )),
     &           a(mpd       ),a(mp(mdl   )),a(mp(mmat  )),
     &           a(mp(mmatside)),a(mp(mmcont)),
     &           a(mp(mdet  )),a(mp(mshl  )),a(mp(mshg  )),
     &           a(mp(mw    )),a(mp(mc    )),a(mpalhs    ),
     &           a(mpdlhs    ),
     &           a(mpbrhs    ),a(mpdiag    ),a(mp(mlm   )),
     &           a(mp(mgrav )),a(mp(mipar )),a(mp(mlado )),
     &           a(mp(mdetc )),a(mp(mshlc )),a(mp(mshgc )),
     &           a(mp(melefd)),a(mp(melred)),a(mp(mshln )),
     &           a(mp(mshgn )),a(mp(mwn   )),a(mp(mdetn )),
     &           a(mp(mdetb )),a(mp(mshlb )),a(mp(mshgb )),
     &           a(mp(mdetpn)),a(mp(mshlpn)),a(mp(mshgpn)),
     &           a(mp(mdside)),a(mp(mxls  )),a(mp(midlsd)),
     &           a(mp(mdsfl )),a(mp(mddis )),a(mp(mdetp )),
     &           a(mp(mshlp )),a(mp(mshgp )),
c     
     &           a(mp(melmbb)),
     &           a(mp(melmcb)),
     &           a(mp(melmcc)),
     &           a(mp(melmbc)),a(mp(melfcc)),
     &           a(mp(melfbb)),a(mp(melfcb)),
     &           a(mp(melmdb)),
c
     &           a(mped      ),

     &           a(mp(mshlpsd)),a(mp(mshlcsd)),
     &           a(mp(mshgpsd)),a(mp(mshgcsd)),
c
     &           a(mp(mshlsd)),
     &           a(mp(mshgsd)),
c
cc     &           a(mp(mepres)),
c     
     &                 numel ,neesq ,nen   ,nsd   ,
     &                 nesd  ,nint  ,neg   ,nrowsh,
     &                 ned   ,nee   ,numnp ,ndof  ,
     &                 ncon  ,nencon,necon ,neep  ,
     &                 nints ,nnods ,nenlad,npars ,
     &                 nside ,nenp  ,nedge ,nodsp ,
     &                 index ,nepcon,nface)
c
      return
c
  400 continue

      call darcy3(a(mp(mien  )),a(mpx       ),a(mp(mxl )),
     &           a(mpd       ),a(mp(mdl   )),a(mp(mmat  )),
     &           a(mp(mmatside)),a(mp(mmcont)),
     &           a(mp(mdet  )),a(mp(mshl  )),a(mp(mshg  )),
     &           a(mp(mw    )),a(mp(mc    )),
     &           a(mp(mgrav )),a(mp(mipar )),a(mp(mlado )),
     &           a(mp(mdetc )),a(mp(mshlc )),a(mp(mshgc )),
     &           a(mp(melefd)),a(mp(melred)),a(mp(mshln )),
     &           a(mp(mshgn )),a(mp(mwn   )),a(mp(mdetn )),
     &           a(mp(mdetb )),a(mp(mshlb )),a(mp(mshgb )),
     &           a(mp(mdetpn)),a(mp(mshlpn)),a(mp(mshgpn)),
     &           a(mp(mdside)),a(mp(mxls  )),a(mp(midlsd)),
     &           a(mp(mdsfl )),a(mp(mddis )),a(mp(mdetp )),
     &           a(mp(mshlp )),a(mp(mshgp )),
c     
     &           a(mp(melmbb)),
     &           a(mp(melmcb)),
     &           a(mp(melmcc)),
     &           a(mp(melfbb)),a(mp(melfcb)),
     &           a(mp(melfcc)),a(mp(melmdb)),
c
     &           a(mped      ),
c
     &           a(mp(mshlpsd)),a(mp(mshlcsd)),
     &           a(mp(mshgpsd)),a(mp(mshgcsd)),
          
     &           a(mp(mshlsd)),
     &           a(mp(mshgsd)),
c
cc     &           a(mp(mepres)),
c
     &           numel ,neesq ,nen   ,
     &           nsd   ,nesd  ,nint  ,
     &           neg   ,nrowsh,ned   ,
     &           nee   ,numnp ,ndof  ,
     &           ncon  ,nencon,necon ,
     &           neep  ,nints ,nnods ,
     &           nenlad, npars,nside ,
     &           nenp  ,nodsp ,index ,
     &           nepcon,nface)
c     
      write(iecho,*)'saiu elastic3'

      if(index.gt.0) then
      call flnorm(a(mp(mien )),a(mpx       ),a(mp(mxl   )),
     &            a(mpd      ),a(mp(mdl   )),a(mp(mmat  )),
     &            a(mp(mmatside)),
     &            a(mp(mc   )),a(mp(mipar )),a(mp(mdlf )) ,
     &            a(mp(mdlp )),a(mp(mdsfl )),a(mp(mdet  )),
     &            a(mp(mshl )),a(mp(mshg  )),a(mp(mw    )),
     &            a(mp(mdetc )),a(mp(mshlc)),a(mp(mshgc )),
     &            a(mp(mddis )),a(mp(mdetp)),a(mp(mshlp )),
     &            a(mp(mshgp )),a(mp(mlado)),
c
     &            a(mped      ),
c
     &            a(mp(mshln )),a(mp(mshgn )),
     &            a(mp(mdetn )),a(mp(mshlb )),a(mp(mshgb )),
     &            a(mp(mdetpn)),a(mp(mshlpn)),a(mp(mshgpn)),
     &            a(mp(mdside)),a(mp(mxls  )),a(mp(midlsd)),
     &            a(mp(mgrav )),a(mp(mwn   )),
c
     &            numel ,neesq ,nen   ,nsd   ,
     &            nesd  ,nint  ,neg   ,nrowsh,
     &            ned   ,nee   ,numnp ,ndof  ,
     &            ncon  ,nencon,necon ,index ,
     &            nints ,iwrite,iout  ,ierrb ,nenp  ,
     &            nside ,nnods ,nenlad,npars ,
     &            nmultp,nodsp)
c
      call flninter(a(mp(mien )),a(mpx       ),a(mp(mxl   )),
     &            a(mpd      ),a(mp(mdl   )),a(mp(mmat  )),
     &            a(mp(mc   )),a(mp(mipar )),a(mp(mdlf )) ,
     &            a(mp(mdlp )),a(mp(mdsfl )),a(mp(mdet  )),
     &            a(mp(mshl )),a(mp(mshg  )),a(mp(mw    )),
     &            a(mp(mdetc )),a(mp(mshlc)),a(mp(mshgc )),
     &            a(mp(mddis )),a(mp(mdetp)),a(mp(mshlp )),
     &            a(mp(mshgp )),
c
     &            a(mp(mshln )),a(mp(mshgn )),
     &            a(mp(mdetn )),a(mp(mshlb )),a(mp(mshgb )),
     &            a(mp(mdetpn)),a(mp(mshlpn)),a(mp(mshgpn)),
     &            a(mp(mdside)),a(mp(mxls  )),a(mp(midlsd)),
     &            a(mp(mgrav )),a(mp(mwn   )),
c
     &            numel ,neesq ,nen   ,nsd   ,
     &            nesd  ,nint  ,neg   ,nrowsh,
     &            ned   ,nee   ,numnp ,ndof  ,
     &            ncon  ,nencon,necon ,index ,
     &            nints ,iwrite,iout1 ,ierrbi,
     &            nenp  ,nside ,nnods ,nenlad,
     &            npars ,nmultp,nodsp )
c
	 end if
c
      return
      end
c***********************************************************************
      subroutine darcy1(shl   ,w   ,       
     &                 c     ,grav  ,
     &                 ien   ,mat   ,
     &                 matside,mcont,
     &                 id    ,lm    ,
     &                 idiag ,ipar  ,
     &                 x     ,lado  ,
     &                 shlc  ,wc    ,
     &                 shlpn ,wpn   ,
     &                 shln  ,wn    ,
     &                 shlb  ,shlp  ,
     &                 wp    ,idside,
c     
     &                 shlpsd,shlcsd,
     &                 shlsd,
c          
     &                 ntype ,numel ,numat ,
     &                 nint  ,nrowsh,nesd  ,
     &                 nen   ,ndof  ,ned   ,
     &                 iprtin,numnp ,ncon  ,
     &                 nencon,necon ,nints ,
     &                 nnods,nnodc  ,nenlad ,npars ,
     &                 nenp ,nside  ,nodsp  ,nedge)

c
c.... program to read, generate and write data for the
c        four-node quadrilateral, elastic continuum element
c
      implicit real*8 (a-h,o-z)
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      dimension shl(nrowsh,nen,*),w(*),
     &          c(9,*),grav(*),ien(nen,*),mat(*),
     &          matside(*),mcont(*),
     &          id(ndof,*),lm(ndof,nodsp,*),idiag(*),
     &          ipar(nodsp,*),idside(nside,*),
     &          x(nesd,*),lado(nside,*),
     &          shlc(nrowsh,nencon,*),wc(*),shlpn(nrowsh,npars,*),
     &          wpn(*),shlp(nrowsh,nenp,*),wp(*),
     &          shln(nrowsh,nnods,*),wn(*),shlb(nrowsh,nenlad,*)
      dimension shlpsd(nrowsh,nenp,*),shlcsd(nrowsh,nencon,*)
      dimension shlsd(nrowsh,nen,*)
c
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
      common /colhtc / neq
c
      write(iecho,1000) ntype,numel,numat
      write(iecho,2000) nint
c
c     geomeria 1D / Lagrange
c
c      generation of local shape functions and weight values
c
c oneshl - pontos de integracao reordenadas: numeracao sequencial(17-11-2012)
c
      call oneshl(shlb,wp,nints,nenlad)
      call oneshl(shln,wn,nints,nnods)
      call oneshl(shlpn,wpn,nints,npars)
c
	if(nen.eq.3) then
      call shlt(shl,w,nint,nen)
      call shlt(shlc,wc,nint,nencon)
      call shlt(shlp,wp,nint,nenp)
c
c shltpbk - pontos de integracao reordendas: numeracao sequencial(17-11-2012)
c
          call shltpbk(shlpsd,nenp,nside,nints)
          call shltpbk(shlcsd,nencon,nside,nints)
          call shltpbk(shlsd,nen,nside,nints)

	else if(nen.eq.4) then
      call shlq(shl,w,nint,nen)
      call shlqpk(shlc,wc,nint,nencon)
      call shlqpk(shlp,wp,nint,nenp)

      call shlqpbk(shlpsd,nenp,nside,nnods,nints)
      call shlqpbk(shlcsd,nencon,nside,nnodc,nints)
      call shlqpbk(shlsd,nen,nside,nenlad,nints)
        else  
           stop
	end if
c
c
      nintb=nside*nints
c
c      read material properties
c
      call fluxmx(c,numat)
c
c	constant body forces
c
      read (iin,5000) (grav(i),i=1,3)
      write (iecho,6000) (grav(i),i=1,3)
c
c    generation of conectivities
c
      call genel(ien,mat,nen)
c
      if (iprtin.eq.0) call prntel(mat,ien,nen,numel)
c
c   generation of conectivety for element multipliers
c
      call genside(idside,nside,nencon)
      call genelad(lado,nside)
      call matlado(lado,mat,matside,mcont,nside,numel,nedge)
      if (iprtin.eq.0) call prntels(mat,lado,nside,numel)
c
      call genelpar(ipar,ien,lado,idside,
     &              nen,nside,nodsp,numel,npars)
c
      if (iprtin.eq.0) call prntelp(mat,ipar,nodsp,numel)
c
c     generation of lm array
c
      call formlm(id,ipar,lm,ndof,ndof,nodsp,numel)
      
      do k=1,numel
         write(iecho,*)'Matriz LM[',k,']'
         do i=1,ndof
           write(iecho,8000)(lm(i,j,k),j=1,nodsp)
         end do
      end do
      write(iecho,8000) 

c
c     modification of idiag array
c
      call colht(idiag,lm,ndof,nodsp,numel,neq)
c      
       write(iecho,*)
     &'IDIAG = mostra a altura das colunas na matriz global'
      write(iecho,9000)(idiag(i),i=1,neq)

c
c     Neuman B.c. not implemented YET
c
c
      return
c
 1000 format(///,
     &' d u a l   h i b r i d   m i x e d   f o r m u l a t i o n ',///
     &//5x,' element type number . . . . . . . . . . .(ntype ) = ',i10
     &//5x,' number of elements  . . . . . . . . . . .(numel ) = ',i10
     &//5x,' number of element material sets . . . . .(numat ) = ',i10)
2000  format(
     &//5x,' numerical integration points  . . . . . .(nint  ) = ',i10)
 5000 format(8f10.0)
 6000 format(////' ',
     &' g r a v i t y   v e c t o r   c o m p o n e n t s      ',//5x,
     &' W-1 direction (sinxsiny)  . . . . . . . = ',      1pe15.8//5x,
     &' W-2 direction (cosxcosy)  . . . . . . . = ',      1pe15.8//5x,
     &' W-3 direction (Polinomial). . . . . . . = ',      1pe15.8//5x)
7000  format(i10)
 8000 format(8(5x,i10))
 9000 format(100(5x,i10))

c
      end
c**** new **********************************************************************
      subroutine fluxmx(c,numat)
c
c.... program to read, write and store properties
c     for plane stres mixed elements
c
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension c(9,*)
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
      do 100 n=1,numat
      if (mod(n,50).eq.1) write(iecho,1000) numat
c
      read(iin,2000) m,del1,del2,del3,del4,del5,del6,del7,del8,del9
      c(1,m)=del1
      c(2,m)=del2
      c(3,m)=del3
      c(4,m)=del4
      c(5,m)=del5
      c(6,m)=del6
      c(7,m)=del7
      c(8,m)=del8
      c(9,m)=del9
      write(iecho,3000) m,del1,del2,del3,del4,del5,del6,del7,del8,del9
c
     
  100 continue
c
      return
c
 1000 format(///,
     &' m a t e r i a l   s e t   d a t a       '   //5x,
     &' number of material sets . . . . . (numat ) = ',i10//,
     & 7x,'set',7x,'del1',7x,'del2',
     & 7x,'del3',7x,'del4', 7x,'del5',7x,'del6',7x,'del7',
     & 7x,'del8',7x,'del9',/)
 2000 format(i10,9f10.0)
 3000 format(i10,2x,9(1x,1pe10.3))
      end
c**** new **********************************************************************
      subroutine darcy0(ien   ,x     ,xl  ,
     &                 d     ,dl    ,mat   ,
     &                 matside,mcont,
     &                 det   ,shl   ,shg   ,
     &                 w     ,c     ,
     &                 grav  ,ipar  ,lado  ,
     &                 detc  ,shlc  ,shgc  ,
     &                 eleffd,elresd,shln  ,
     &                 shgn  ,wn    ,detn  ,
     &                 detb  ,shlb  ,shgb  ,
     &                 detpn ,shlpn ,shgpn ,
     &                 idside,xls   ,idlsd ,
     &                 dsfl  ,ddis  ,detp  ,
     &                 shlp  ,shgp  ,
     &                 elmbb ,
     &                 elmcb ,
     &                 elfbb ,elfcb ,
     &                 elmdb ,
c     
     &                 ideg  ,
c
     &                 shlpsd,shlcsd,
     &                 shgpsd,shgcsd,
c     
     &                 shlsd,
     &                 shgsd,
c     
     &                 numel ,neesq ,nen   ,
     &                 nsd   ,nesd  ,nint  ,
     &                 neg   ,nrowsh,ned   ,
     &                 nee   ,numnp ,ndof  ,
     &                 ncon  ,nencon,necon ,
     &                 neep  ,nints ,nnods ,
     &                 nenlad,npars ,nside ,
     &                 nenp  ,nodsp ,index ,
     &                 nepcon,nface )
c
c
c.... program to calculate stifness matrix and force array for the
c        plane elasticity element and
c        assemble into the global left-hand-side matrix
c        and right-hand side vector
c
      implicit real*8 (a-h,o-z)
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      logical diag,quad,zerodl
cc      dimension aut1(10),aut2(10),un1(10),un2(10)
      dimension elmbb(nepcon,*),elmcb(nepcon,*)
      dimension elmdb(nee,*)
      dimension elfbb(*),elfcb(*)
      dimension ien(nen,*),x(nsd,*),xl(nesd,*),d(ndof,*),dl(ndof,*),
     &          mat(*),det(*),shl(nrowsh,nen,*),shg(nrowsh,nen,*),
     &          w(*),c(9,*),matside(*),mcont(*),
     &          grav(*),ipar(nodsp,*),lado(nside,*),detc(*),
     &          shlc(nrowsh,nencon,*),shgc(nrowsh,nencon,*),
     &          eleffd(nee,*),elresd(*)
      dimension shln(3,nnods,*),wn(*),detn(*),shgn(3,nnods,*),
     &          detb(*),shlb(3,nenlad,*),shgb(3,nenlad,*),
     &          detpn(*),shlpn(3,npars,*),shgpn(3,npars,*),
     &          idside(nside,*),xls(nesd,*),idlsd(*),
     &          dsfl(ncon,nencon,*),ddis(ned,nenp,*) 
      dimension dls(ndof,12),detp(*),shlp(3,nenp,*),shgp(3,nenp,*)
      dimension ue(6),duex(6),duey(6)
c
      dimension ideg(2*ndof,*)
c
      dimension shlpsd(nrowsh,nenp,*),shlcsd(nrowsh,nencon,*),
     &          shgpsd(nrowsh,nenp,*),shgcsd(nrowsh,nencon,*)
c
      dimension shlsd(nrowsh,nen,*),
     &          shgsd(nrowsh,nen,*)
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
c      consistent matrix
c
      diag = .false.
      pi=4.d00*datan(1.d00)
      gf1=grav(1)
      gf2=grav(2)
      gf3=grav(3)      
c
c
      do 500 nel=1,numel
c
c      clear stiffness matrix and force array
c
      call clear(elmbb,nepcon*nepcon)
      call clear(elfbb,nepcon)
c
c      localize coordinates and Dirichlet b.c.
c
      call local(ien(1,nel),x,xl,nen,nsd,nesd)
c
c
      m = mat(nel)
      quad = .true.
      if (nen.eq.4.and.ien(3,nel).eq.ien(4,nel)) quad = .false.
c     
c
      nintb=nside*nints
      call shgqs(xl,detc,shlc,shgc,nint,nel,neg,quad,nencon,shl,nen)
      call shgqs(xl,detp,shlp,shgp,nint,nel,neg,quad,nenp,shl,nen)
c
      call shgtqsd(ien(1,nel),xl,shlpsd,shgpsd,
     &     nside,nintb,nints,nel,neg,nenp,shlsd,nen)
      call shgtqsd(ien(1,nel),xl,shlcsd,shgcsd,
     &     nside,nintb,nints,nel,neg,nencon,shlsd,nen)
c
c....... form stiffness matrix
c
c      calculo de h - caracteristico para cada elemento
c
      if(nen.eq.3) then
c
        h2=xl(1,2)*xl(2,3)+xl(1,1)*xl(2,2)+xl(1,3)*xl(2,1)
     &    -xl(1,1)*xl(2,3)-xl(1,2)*xl(2,1)-xl(1,3)*xl(2,2)
        h2=h2*pt5
        h=dsqrt(h2)
c
      else if (nen.eq.4) then
c
        h2=0
        do 100 l=1,2
         h2=h2+(xl(1,l)-xl(1,l+2))**2+(xl(2,l)-xl(2,l+2))**2
 100    continue 
        h=dsqrt(h2)/2.d00
        h2=h*h
c
      else 
        stop
      end if
c
c      set up material properties
c     
	  eps=c(3,m)  !
 	  delta1=c(4,m)
 	  delta2=c(5,m)
	  del1=c(6,m)
 	  del2=c(7,m)
        del3=c(8,m)
c
c     Matriz de condutividade anisotropica
c
c            | xka1  xka2 |
c K1 = alpha*|            | ,
c            | xka2  xka1 |
c
        alpha=c(9,m)
	  xka1=c(1,m)*alpha
	  xka2=c(2,m)*alpha
c
        xkl = c(1,m)*alpha
c
        betap = xkl*delta1*h**delta2     !betap
c
c
c...  Matriz de Resistencia Hidraulica ({K}-1)
c
c               | a11  a12 |
c inv{K1} = A = |          | ,
c               | a21  a22 |
c
      deta=xka1*xka1-xka2*xka2
      a11 = xka1/deta
      a12 =-xka2/deta
      a22 = a11
      a21 = a12
c
c
        if(nel.eq.1) then
         write(iecho,*)'xka1,xka2,eps'
         write(iecho,*)xka1,xka2,eps
         write(iecho,*)'delta1,delta2,betah'
         write(iecho,*)delta1,delta2,betah
         write(iecho,*)'delta3,delta4,betan'
         write(iecho,*)delta3,delta4,betan
        endif
c
c
c.....loop on integration points
c
      do 400 l=1,nint
      c1=detp(l)*w(l)
c      
      xx=0.d00
      yy=0.d00
c
      do 2323 i=1,nen
      xx=xx+shl(3,i,l)*xl(1,i)
      yy=yy+shl(3,i,l)*xl(2,i)
2323  continue
c
      pix=pi*xx
      piy=pi*yy
      sx=dsin(pix)
      sy=dsin(piy)
      cx=dcos(pix)
      cy=dcos(piy)
c
      pi2=pi*pi
c
c
c     fonte Darcy
c
      if (index.eq.1) then
c
        if(m.eq.1) then
            fd = 4.d00*pi2*sx*sy
c
        else
            fd = 2.d00*pi2*alpha*(2.d00*sx*sy - cx*cy)
c
        end if
c
      else if(index.eq.2) then
c
      sx=dsin(xx)
      sy=dsin(yy)
      cx=dcos(xx)
      cy=dcos(yy)
      ex=dexp(xx)

        if(m.eq.1) then
            fd = alpha*xx*(2.d00*sy+cy) + sy
c
        else
            fd = -2.d00*alpha*ex*cy
c
        end if
c
      else if(index.eq.3) then
            fd = cx*cy
c
      end if
c
c      loop to compute volume integrals
c
      do 3000 j=1,nencon
      djx=shgc(1,j,l)*c1
      djy=shgc(2,j,l)*c1
      djn=shgc(3,j,l)*c1
c
c     curl(Av)
c
      rtj1=a11*djy - a21*djx
      rtj2=a12*djy - a22*djx
c
c.... source terms      
c
      naj=ncon*(j-1)
      naj1=naj+1
      naj2=naj+2
c
c     (f,div v)/k1
c
c   v1*f1
      elfbb(naj1)=elfbb(naj1) + del1*djx*fd/xkl
c   v2*f2
      elfbb(naj2)=elfbb(naj2) + del1*djy*fd/xkl
c
c
c.... element stiffness
c
      do 3003 i=1,nencon
      nai=ncon*(i-1)
      nai1=nai+1
      nai2=nai+2
c
      dix=shgc(1,i,l)
      diy=shgc(2,i,l)
      din=shgc(3,i,l)
c
c     curl(Au)
c
      rti1=a11*diy - a21*dix
      rti2=a12*diy - a22*dix
c
c
c termos:  (Av,u)
c
c
      elmbb(naj1,nai1) = elmbb(naj1,nai1)
     &                 + a11*djn*din
c
      elmbb(naj1,nai2) = elmbb(naj1,nai2)
     &                 + a12*djn*din
c
      elmbb(naj2,nai1) = elmbb(naj2,nai1)
     &                 + a21*djn*din
c
      elmbb(naj2,nai2) = elmbb(naj2,nai2) 
     &                 + a22*djn*din
c
c termos:  del3*K(Av,Au)
c
c
      elmbb(naj1,nai1) = elmbb(naj1,nai1)
     &                 + del3*djn*din*a11
c
      elmbb(naj1,nai2) = elmbb(naj1,nai2)
     &                 + del3*djn*din*a12
c
      elmbb(naj2,nai1) = elmbb(naj2,nai1)
     &                 + del3*djn*din*a21
c
      elmbb(naj2,nai2) = elmbb(naj2,nai2)
     &                 + del3*djn*din*a22
c
c
c termos:  del1*(div(v), div(u))/k1
c
c     v1*u1: 
      elmbb(naj1,nai1) = elmbb(naj1,nai1) 
     &                 + del1*djx*dix/xkl
c     v2*u1:
      elmbb(naj2,nai1) = elmbb(naj2,nai1)
     &                 + del1*djy*dix/xkl
c     v1*u2:
      elmbb(naj1,nai2) = elmbb(naj1,nai2)
     &                 + del1*djx*diy/xkl
c     v2*u2:
      elmbb(naj2,nai2) = elmbb(naj2,nai2) 
     &                 + del1*djy*diy/xkl
c
c	  
c
c termos: del2*(k1 rot(Av), rot(Au))
c
c     v1*u1: 
      elmbb(naj1,nai1) = elmbb(naj1,nai1) 
     &                 + del2*rti1*rtj1*xkl
c     v2*u1:
      elmbb(naj2,nai1) = elmbb(naj2,nai1)
     &                 + del2*rtj2*rti1*xkl
c     v1*u2:
      elmbb(naj1,nai2) = elmbb(naj1,nai2)
     &                 + del2*rtj1*rti2*xkl
c     v2*u2:
      elmbb(naj2,nai2) = elmbb(naj2,nai2) 
     &                 + del2*rti2*rtj2*xkl
c
 3003 continue
c
      do 3005 i=1,nenp
	  nbi1=i+necon
c       
      dix=shgp(1,i,l)
      diy=shgp(2,i,l)
      din=shgp(3,i,l)
c
c
c termos:  -(p, div(v)) 
c
c     p*v1: 
      elmbb(naj1,nbi1) = elmbb(naj1,nbi1) - din*djx 
c
c     p*v2: 
      elmbb(naj2,nbi1) = elmbb(naj2,nbi1) - din*djy
c
c
c termos:  del3*K*(gradp, Av)
c
c     p*v1: 
      elmbb(naj1,nbi1) = elmbb(naj1,nbi1) + del3*dix*djn
c
c     p*v2: 
      elmbb(naj2,nbi1) = elmbb(naj2,nbi1) + del3*diy*djn
c
c
 3005 continue
 3000 continue
c
c
      do 303 j=1,nenp
	  nbj1=j+necon
      djx=shgp(1,j,l)*c1
      djy=shgp(2,j,l)*c1
      djn=shgp(3,j,l)*c1
c
c     - (f,q)
c
      elfbb(nbj1)=elfbb(nbj1) - fd*djn
c
c
      do 304 i=1,nencon
      nai=ncon*(i-1)
      nai1=nai+1
      nai2=nai+2
c
      dix=shgc(1,i,l)
      diy=shgc(2,i,l)
      din=shgc(3,i,l)
c
c
c termos:  -(div(u),q) 
c
c     u1*q: 
      elmbb(nbj1,nai1) = elmbb(nbj1,nai1) - djn*dix 
c
c
c     u2*q: 
      elmbb(nbj1,nai2) = elmbb(nbj1,nai2) - djn*diy
c
c
c termos:  del3*K*(Au,gradq)
c
c     u1*q: 
      elmbb(nbj1,nai1) = elmbb(nbj1,nai1) + del3*djx*din
c
c
c     u2*q: 
      elmbb(nbj1,nai2) = elmbb(nbj1,nai2) + del3*djy*din
c
c
  304 continue
      do 305 i=1,nenp
      nbi1=i+necon     
      dix=shgp(1,i,l)
      diy=shgp(2,i,l)
      din=shgp(3,i,l)
c
c     termo: - eps*(p,q)
c
c     p*q: 
      elmbb(nbj1,nbi1) = elmbb(nbj1,nbi1) - eps*din*djn
c
c
c
c termos:  del3*K*(gradp,gradq)
c
      elmbb(nbj1,nbi1) = elmbb(nbj1,nbi1)
     &                 + del3*xka1*(dix*djx + diy*djy)
     &                 + del3*xka2*(diy*djx + dix*djy)
c
c
  305 continue
  303 continue
c
  400 continue
c
c
c......Symmetrization - boundary terms 
c
      do 4000 ns=1,nside
c
c-----localiza os no's do lado ns
c
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,nel)
      nl2=ien(ns2,nel)
c
      if(nl2.gt.nl1) then 
	  sign = 1.d00
        do  nn=1,nenlad
           idlsd(nn)=idside(ns,nn)
	  end do
      else
	  sign = -1.d00
	  idlsd(1) = idside(ns,2)
	  idlsd(2) = idside(ns,1)
	  id3 = nenlad-2
	  if(id3.gt.0) then
	    do il=3,nenlad
	      idlsd(il) = idside(ns,nenlad+3-il)
	    end do 
	  end if
      end if 
c
c    localize the coordinates of the side nodes 
c
      do 2000 nn=1,nenlad
      nl=idlsd(nn)
      xls(1,nn)=xl(1,nl)
      xls(2,nn)=xl(2,nl)
 2000 continue
c
      call oneshgp(xls,detn,shlb,shln,shgn,
     &             nenlad,nnods,nints,nesd,ns,nel,neg) 
      call oneshgp(xls,detpn,shlb,shlpn,shgpn,
     &             nenlad,npars,nints,nesd,ns,nel,neg) 
c
c     projecao local do multiplicador
c
      ndgs = lado(ns,nel) ! numero global do lado
c
cc      do igl = 1
        igl = 1
        call drchbc(shgpn,shlb,detpn,wn,
     &            xls,dls,pi,eps,alpha,ndgs,
     &            nints,nenlad,npars,ndof,igl,m,index)
	ngs = npars*(ndgs-1)
	nls = npars*(ns-1)
        do i=1,npars
          d(igl,ngs + i) = dls(igl,i)
        end do
c
c.....compute boundary integral - symmetrization
c
      do 1000 ls = 1,nints
      lb = (ns-1)*nints + ls
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c
c    valor exato do multiplicador
c

         call uexafx(x1,x2,ue,duex,duey,alphas,m,index)
c
c
       dhs=ue(3)
c
c
       do 1110 j=1,nencon
          naj = ncon*(j-1)
          naj1=naj+1
          naj2=naj+2
          djn=shgcsd(3,j,lb)*detn(ls)*wn(ls)
          djx=shgcsd(1,j,lb)*detn(ls)*wn(ls)
          djy=shgcsd(2,j,lb)*detn(ls)*wn(ls)
          gjn=djx*xn1 + djy*xn2
c
c         termo: -(p_exato,v.n)
c
          elfbb(naj1)= elfbb(naj1) - dhs*djn*xn1
c
          elfbb(naj2)= elfbb(naj2) - dhs*djn*xn2
c
c
 1110 continue
c
c    jump terms
c
       do 9100 j=1,nenp
          nbj1=j+necon
          djn=shgpsd(3,j,lb)*detn(ls)*wn(ls)
          djx=shgpsd(1,j,lb)*detn(ls)*wn(ls)
          djy=shgpsd(2,j,lb)*detn(ls)*wn(ls)
          gjn=djx*xn1 + djy*xn2
          djn1=djn*xn1
          djn2=djn*xn2
c
c         termo:  betap*(q,p_exato)
c
          elfbb(nbj1) = elfbb(nbj1) + betap*dhs*djn
c
       do 9200 i=1,nenp
          nbi1=i+necon
          din=shgpsd(3,i,lb)
          dix=shgpsd(1,i,lb)
          diy=shgpsd(2,i,lb)
          gin=dix*xn1 + diy*xn2
          din1=din*xn1
          din2=din*xn2
c
c         termo:  betap*(p,q)
c
          elmbb(nbj1,nbi1) = elmbb(nbj1,nbi1) + betap*din*djn
c 
 9200 continue
 9100 continue
c
 1000 continue
 4000 continue
c
c	Local TLDG-solution (potencial) 
c
c
      if(nel.eq.1) then  
c	 write(25,*) nel,nepcon
       do i=1,nepcon
c         write(25,91) (elmbb(i,j),j=1,nepcon),elfbb(i)
c        write(25,91) (elmbb(j,j),j=1,nepcon)
       end do 
       end if
  91   format(30e15.5)
c
      call solvetdg(elmbb,elfbb,nepcon)
c          
c          
c      valores nodais descontinuos
c
        jk=0
        do 9171 j=1,nencon
        do 9272 jj=1,ncon
        jk=jk+1
        dsfl(jj,j,nel)=elfbb(jk)
 9272   continue
 9171   continue
c
      do j=1,nenp
	  jn=j+necon
	  ddis(1,j,nel) = elfbb(jn)
      end do
c
  500 continue
c
      return
      end
      
c**** new **********************************************************************
      subroutine darcy2(ien   ,x     ,xl  ,
     &                 d     ,dl    ,mat   ,
     &                 matside,mcont,
     &                 det   ,shl   ,shg   ,
     &                 w     ,c     ,alhs  ,
     &                 dlhs  ,
     &                 brhs  ,idiag ,lm    ,
     &                 grav  ,ipar  ,lado  ,
     &                 detc  ,shlc  ,shgc  ,
     &                 eleffd,elresd,shln  ,
     &                 shgn  ,wn    ,detn  ,
     &                 detb  ,shlb  ,shgb  ,
     &                 detpn ,shlpn ,shgpn ,
     &                 idside,xls   ,idlsd ,
     &                 dsfl  ,ddis  ,detp  ,
     &                 shlp  ,shgp  ,
     &                 elmbb ,
     &                 elmcb ,elmcc ,
     &                 elmbc ,elfcc ,
     &                 elfbb ,elfcb ,
     &                 elmdb ,
c     
     &                 ideg  ,
c     
     &                 shlpsd,shlcsd,
     &                 shgpsd,shgcsd,
c     
     &                 shlsd,
     &                 shgsd,
c
cc     &                 epres,
c     
     &                 numel ,neesq ,nen   ,
     &                 nsd   ,nesd  ,nint  ,
     &                 neg   ,nrowsh,ned   ,
     &                 nee   ,numnp ,ndof  ,
     &                 ncon  ,nencon,necon ,
     &                 neep  ,nints ,nnods ,
     &                 nenlad,npars ,nside ,
     &                 nenp  ,nedge ,nodsp ,
     &                 index ,nepcon,nface)
c
c
c.... program to calculate stifness matrix and force array for the
c        plane elasticity element and
c        assemble into the global left-hand-side matrix
c        and right-hand side vector
c
      implicit real*8 (a-h,o-z)
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      logical diag,quad,zerodl
      dimension elmbb(nepcon,*),elmcb(nepcon,*),elmbc(nee,*)
      dimension elmdb(nee,*),elmcc(nee,*)
      dimension elfbb(*),elfcb(*),elfcc(*)
      dimension ien(nen,*),x(nsd,*),xl(nesd,*),d(ndof,*),dl(ndof,*),
     &          mat(*),det(*),shl(nrowsh,nen,*),shg(nrowsh,nen,*),
     &          w(*),c(9,*),alhs(*),brhs(*),idiag(*),lm(ndof,nodsp,*),
     &          grav(*),ipar(nodsp,*),lado(nside,*),matside(*),
     &          detc(*),shlc(nrowsh,nencon,*),shgc(nrowsh,nencon,*),
     &          eleffd(nee,*),elresd(*),dlhs(*),mcont(*)
      dimension shln(3,nnods,*),wn(*),detn(*),shgn(3,nnods,*),
     &          detb(*),shlb(3,nenlad,*),shgb(3,nenlad,*),
     &          detpn(*),shlpn(3,npars,*),shgpn(3,npars,*),
     &          idside(nside,*),xls(nesd,*),idlsd(*),
     &          dsfl(ncon,nencon,*),ddis(ned,nenp,*) 
	  dimension dls(ndof,12),detp(*),shlp(3,nenp,*),shgp(3,nenp,*)
c
      dimension ideg(2*ndof,*)

      dimension shlpsd(nrowsh,nenp,*),shlcsd(nrowsh,nencon,*),
     &          shgpsd(nrowsh,nenp,*),shgcsd(nrowsh,nencon,*)
c
      dimension shlsd(nrowsh,nen,*),
     &          shgsd(nrowsh,nen,*)
      dimension imprime(nee,nee)
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
c      consistent matrix
c
      diag = .false.
      pi=4.d00*datan(1.d00)
      gf1=grav(1)
      gf2=grav(2)
      gf3=grav(3)      
c
      
      do 500 nel=1,numel
c
c      clear stiffness matrix and force array
c
      call clear(eleffd,nee*nee)
      call clear(elmbb,nepcon*nepcon)
      call clear(elmcb,nepcon*nee)
      call clear(elmbc,nee*nepcon)
      call clear(elmdb,nee*nepcon)
c      
      call clear(elfbb,nepcon)
	  call clear(elfcb,nee)
	  call clear(elresd,nee)

c
c     Multiplier matrix
c
      call clear(eleffd,nee*nee)
      call clear(elresd,nee)
c
c      localize coordinates and Dirichlet b.c.
c
      call local(ien(1,nel),x,xl,nen,nsd,nesd)
      call local(ipar(1,nel),d,dl,nodsp,ndof,ndof)
c
c
      m = mat(nel)
      quad = .true.
      if (nen.eq.4.and.ien(3,nel).eq.ien(4,nel)) quad = .false.
c     
c
      call shgqs(xl,detc,shlc,shgc,nint,nel,neg,quad,nencon,shl,nen)
      call shgqs(xl,detp,shlp,shgp,nint,nel,neg,quad,nenp,shl,nen)
c
      nintb=nside*nints
c      
      call shgtqsd(ien(1,nel),xl,shlpsd,shgpsd,
     &     nside,nintb,nints,nel,neg,nenp,shlsd,nen)
      call shgtqsd(ien(1,nel),xl,shlcsd,shgcsd,
     &     nside,nintb,nints,nel,neg,nencon,shlsd,nen)

c
c....... form stiffness matrix
c
c      calculo de h - caracteristico para cada elemento
c
      if(nen.eq.3) then
c
        h2=xl(1,2)*xl(2,3)+xl(1,1)*xl(2,2)+xl(1,3)*xl(2,1)
     &    -xl(1,1)*xl(2,3)-xl(1,2)*xl(2,1)-xl(1,3)*xl(2,2)
        h2=h2*pt5
        h=dsqrt(h2)
c
      else if (nen.eq.4) then
c
        h2=0
        do 100 l=1,2
         h2=h2+(xl(1,l)-xl(1,l+2))**2+(xl(2,l)-xl(2,l+2))**2
 100    continue 
        h=dsqrt(h2)/2.d00
        h2=h*h
c
      else 
        stop
      end if

c------------------------fonte Darcy tracador---------------------------
c
      if (nel.eq.1) then
        fd = 200.0d00
        
      else if (nel.eq.numel) then
        fd = -200.0d00
        
      else
        fd = 0.0d00
        
      end if
c-----------------------------------------------------------------------


c      set up material properties
c     
	  eps=c(3,m)  !
 	  delta1=c(4,m)
 	  delta2=c(5,m)
	  del1=c(6,m)
 	  del2=c(7,m)
        del3=c(8,m)
c
c     Matriz de condutividade anisotropica
c
c            | xka1  xka2 |
c K1 = alpha*|            | ,
c            | xka2  xka1 |
c



c----------interpolacao concentracao tracador----------
cca=0.d0
c      do 201 j=1,nenc
c      cca= cca+ shgc(3,j,l)*ccl(1,j)
c201   continue
c
c     Para garantir que a concentracao fique sempre entre 0 e 1
c
c      if (cca.lt.0.d00) then
!       write(60,*) cca, nustep
c      cca=0.d00
c      endif 
      
c      if (cca.gt.1.d00) then
!       write(60,*) ccca
c      cca=1.d00
c      endif
      
c      rmu =rmu0/((1.d0 + cca*((rmm**0.25d0)-1.d0))**4)




        alpha=c(9,m)
	  xka1=c(1,m)*alpha
	  xka2=c(2,m)*alpha
c
c
c...  Matriz de Resistencia Hidraulica ({K}-1)
c
c               | a11  a12 |
c inv{K1} = A = |          | ,
c               | a21  a22 |
c
      deta=xka1*xka1-xka2*xka2
      a11 = xka1/deta
      a12 =-xka2/deta
      a22 = a11
      a21 = a12
c
c
      xkl = c(1,m)*alpha
c
        betap = xkl*delta1*h**delta2      !betap
c
c.....loop on integration points
c
      do 400 l=1,nint
      c1=detp(l)*w(l)

c      
      xx=0.d00
      yy=0.d00
c
      do 2323 i=1,nen
      xx=xx+shl(3,i,l)*xl(1,i)
      yy=yy+shl(3,i,l)*xl(2,i)
2323  continue
c
      pix=pi*xx
      piy=pi*yy
      sx=dsin(pix)
      sy=dsin(piy)
      cx=dcos(pix)
      cy=dcos(piy)
c
      pi2=pi*pi
c
c
c     fonte Darcy 
c-----------------------------------------------------------------------
c      if (index.eq.1) then
c
c        if(m.eq.1) then
c            fd = 4.d00*pi2*sx*sy
c
c        else
c            fd = 2.d00*pi2*alpha*(2.d00*sx*sy - cx*cy)
c
c        end if
c
c      else if(index.eq.2) then
c
c      sx=dsin(xx)
c      sy=dsin(yy)
c      cx=dcos(xx)
c      cy=dcos(yy)
c      ex=dexp(xx)
c
c        if(m.eq.1) then
c            fd = alpha*xx*(2.d00*sy+cy) + sy
c
c        else
c            fd = -2.d00*alpha*ex*cy
c
c        end if
c
c      else if(index.eq.3) then
c            fd = cx*cy
c
c      end if
c-----------------------------------------------------------------------

c
c      loop to compute volume integrals
c
      do 3000 j=1,nencon
      djx=shgc(1,j,l)*c1
      djy=shgc(2,j,l)*c1
      djn=shgc(3,j,l)*c1
c
c     curl(Av)
c
      rtj1=a11*djy - a21*djx
      rtj2=a12*djy - a22*djx
c
c.... source terms      
c
      naj=ncon*(j-1)
      naj1=naj+1
      naj2=naj+2
c
c     (f,div v)/k1
c
      elfbb(naj1)=elfbb(naj1) + del1*djx*fd/xkl
      elfbb(naj2)=elfbb(naj2) + del1*djy*fd/xkl
c
c
c
c.... element stiffness
c
      do 3003 i=1,nencon
      nai=ncon*(i-1)
      nai1=nai+1
      nai2=nai+2
c
      dix=shgc(1,i,l)
      diy=shgc(2,i,l)
      din=shgc(3,i,l)
c
c     curl(Au)
c
      rti1=a11*diy - a21*dix
      rti2=a12*diy - a22*dix
c
c
c termos:  (Av,u)
c
c
      elmbb(naj1,nai1) = elmbb(naj1,nai1) + a11*djn*din
      elmbb(naj1,nai2) = elmbb(naj1,nai2) + a12*djn*din
c
      elmbb(naj2,nai1) = elmbb(naj2,nai1) + a21*djn*din
      elmbb(naj2,nai2) = elmbb(naj2,nai2) + a22*djn*din
c
c termos:  del3*K(Av,Au)
c
c
      elmbb(naj1,nai1) = elmbb(naj1,nai1) + del3*djn*din*a11
      elmbb(naj1,nai2) = elmbb(naj1,nai2) + del3*djn*din*a12
c
      elmbb(naj2,nai1) = elmbb(naj2,nai1) + del3*djn*din*a21
      elmbb(naj2,nai2) = elmbb(naj2,nai2) + del3*djn*din*a22
c
c
c termos:  del1*(div(v), div(u))/k1
c
      elmbb(naj1,nai1) = elmbb(naj1,nai1) + del1*djx*dix/xkl
      elmbb(naj1,nai2) = elmbb(naj1,nai2) + del1*djx*diy/xkl
c
      elmbb(naj2,nai1) = elmbb(naj2,nai1) + del1*djy*dix/xkl
      elmbb(naj2,nai2) = elmbb(naj2,nai2) + del1*djy*diy/xkl
c
c	  
c
c termos: del2*(k1 rot(Av), rot(Au))
c
      elmbb(naj1,nai1) = elmbb(naj1,nai1) + del2*rtj1*rti1*xkl
      elmbb(naj1,nai2) = elmbb(naj1,nai2) + del2*rtj1*rti2*xkl
c     
      elmbb(naj2,nai1) = elmbb(naj2,nai1) + del2*rtj2*rti1*xkl
      elmbb(naj2,nai2) = elmbb(naj2,nai2) + del2*rtj2*rti2*xkl
c
 3003 continue
c
      do 3005 i=1,nenp
	  nbi1=i+necon
c       
      dix=shgp(1,i,l)
      diy=shgp(2,i,l)
      din=shgp(3,i,l)
c
c
c termos:  -(p, div(v)) 
c
      elmbb(naj1,nbi1) = elmbb(naj1,nbi1) - din*djx 
      elmbb(naj2,nbi1) = elmbb(naj2,nbi1) - din*djy
c
c
c termos:  del3*K*(gradp, Av)
c
      elmbb(naj1,nbi1) = elmbb(naj1,nbi1) + del3*dix*djn
      elmbb(naj2,nbi1) = elmbb(naj2,nbi1) + del3*diy*djn
c
c
 3005 continue
 3000 continue
c
c
c
      do 303 j=1,nenp
	  nbj1=j+necon
      djx=shgp(1,j,l)*c1
      djy=shgp(2,j,l)*c1
      djn=shgp(3,j,l)*c1
c
c     (f,q)
c
      elfbb(nbj1)=elfbb(nbj1) - fd*djn
c
c
      do 304 i=1,nencon
      nai=ncon*(i-1)
      nai1=nai+1
      nai2=nai+2
c
      dix=shgc(1,i,l)
      diy=shgc(2,i,l)
      din=shgc(3,i,l)
c
c
c termos:  -(div(u),q) 
c
      elmbb(nbj1,nai1) = elmbb(nbj1,nai1) - djn*dix 
      elmbb(nbj1,nai2) = elmbb(nbj1,nai2) - djn*diy
c
c
c termos:  del3*K*(Au,gradq)
c
      elmbb(nbj1,nai1) = elmbb(nbj1,nai1) + del3*djx*din
      elmbb(nbj1,nai2) = elmbb(nbj1,nai2) + del3*djy*din
c
c
  304 continue
      do 305 i=1,nenp
	  nbi1=i+necon     
      dix=shgp(1,i,l)
      diy=shgp(2,i,l)
      din=shgp(3,i,l)
c
c     termo: - eps*(p,q)
c
      elmbb(nbj1,nbi1) = elmbb(nbj1,nbi1) - eps*din*djn
c
c
c
c termos:  del3*K*(gradp,gradq)
c
      elmbb(nbj1,nbi1) = elmbb(nbj1,nbi1)
     &                 + del3*xka1*(dix*djx + diy*djy)
     &                 + del3*xka2*(diy*djx + dix*djy)
c
c
  305 continue
  303 continue
c
  400 continue
c
c......boundary terms - due to symmetrization
c      typical of disconyinuous Galerkin methods 
c      with SIP 
c
      do 4000 ns=1,nside
c
c-----localiza os no's do lado ns
c
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,nel)
      nl2=ien(ns2,nel)
c      
      if(nl2.gt.nl1) then 
	  sign = 1.d00
        do  nn=1,nenlad
           idlsd(nn)=idside(ns,nn)
	  end do
	else
	  sign = -1.d00
	  idlsd(1) = idside(ns,2)
	  idlsd(2) = idside(ns,1)
	  id3 = nenlad-2
	  if(id3.gt.0) then
	    do il=3,nenlad
	      idlsd(il) = idside(ns,nenlad+3-il)
	    end do 
	  end if
      end if 
c
c    localize the coordinates of the side nodes 
c
      do 2000 nn=1,nenlad
      nl=idlsd(nn)
      xls(1,nn)=xl(1,nl)
      xls(2,nn)=xl(2,nl)
 2000 continue
c
      call oneshgp(xls,detn,shlb,shln,shgn,
     &             nenlad,nnods,nints,nesd,ns,nel,neg)     
      call oneshgp(xls,detpn,shlb,shlpn,shgpn,
     &             nenlad,npars,nints,nesd,ns,nel,neg) 
c
c     Dirichlet boundary condition
c
      ndgs = lado(ns,nel) ! numero global do lado
      igl = 1
        if(ideg(igl,ndgs).eq.1) then
        call drchbc(shgpn,shlb,detpn,wn,
     &            xls,dls,pi,eps,alpha,ndgs,
     &            nints,nenlad,npars,ndof,igl,m,index)
	ngs = npars*(ndgs-1)
	nls = npars*(ns-1)
        do i=1,npars
          d(igl,ngs + i) = dls(igl,i)
	    dl(igl,nls + i) = dls(igl,i)
        end do
        end if
c
c     Neumann boundary condition (lambda_p)
c
	if(ideg(2,ndgs).eq.1) then
      call neuhbc(shgpn,shlb,detpn,wn,idlsd,
     &            xls,elresd,pi,eps,sign,ndof,alpha,m,
     &            nints,nenlad,ns,npars,index)
	end if
c
c.....compute boundary integral - symmetrization
c
      do 1000 ls = 1,nints
c
      lb = (ns-1)*nints + ls
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
c      Calculo do vetor normal
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c
       do 1100 j=1,nencon
          naj = ncon*(j-1)
          naj1=naj+1
	      naj2=naj+2
          djn=shgcsd(3,j,lb)*detn(ls)*wn(ls)
          djx=shgcsd(1,j,lb)*detn(ls)*wn(ls)
          djy=shgcsd(2,j,lb)*detn(ls)*wn(ls)
c
c
	   do i=1,npars
          nci = (ns-1)*npars  + i
		  din = shgpn(2,i,ls)
c
c         termo: + (v.n,\lambda_p)
c
          elmcb(naj1,nci)= elmcb(naj1,nci)
     &                    + din*xn1*djn 
c
          elmcb(naj2,nci)= elmcb(naj2,nci)
     &                    + din*xn2*djn 	
c
c         termo: + (u.n,\mu_p)
c
          elmbc(nci,naj1)= elmbc(nci,naj1)
     &                    + din*xn1*djn 
c
          elmbc(nci,naj2)= elmbc(nci,naj2)
     &                    + din*xn2*djn 	
         end do
c
 1100 continue
c
c      jump terms
c
c 
       do 9100 j=1,nenp
          nbj1=j+necon
          djn=shgpsd(3,j,lb)*detn(ls)*wn(ls)
          djx=shgpsd(1,j,lb)*detn(ls)*wn(ls)
          djy=shgpsd(2,j,lb)*detn(ls)*wn(ls)
c
       do 9200 i=1,nenp
          nbi1=i+necon
          din=shgpsd(3,i,lb)
          dix=shgpsd(1,i,lb)
          diy=shgpsd(2,i,lb)
c
c         termo:  betap*(p,q)
c
          elmbb(nbj1,nbi1) = elmbb(nbj1,nbi1) + betap*din*djn
c 
 9200 continue
 9100 continue
c
c
       do 8100 j=1,nenp
          nbj1=j+necon
          djn=shgpsd(3,j,lb)*detn(ls)*wn(ls)
          djx=shgpsd(1,j,lb)*detn(ls)*wn(ls)
          djy=shgpsd(2,j,lb)*detn(ls)*wn(ls)
c
          do i=1,npars
            nci = (ns-1)*npars  + i
            din = shgpn(2,i,ls)
c
c         termo: - betap*(q,\lambda_p)
c
          elmcb(nbj1,nci)= elmcb(nbj1,nci) - betap*din*djn
c
c
c         termo: - betap*(p,\mu_p)
c
          elmbc(nci,nbj1)= elmbc(nci,nbj1) - betap*din*djn
c          
       end do
 8100 continue
c
        do 8111 j=1,npars
          ncj = (ns-1)*npars  + j
          djn=shgpn(2,j,ls)*detpn(ls)*wn(ls)
c
          do i=1,npars
            nci = (ns-1)*npars  + i
            din = shgpn(2,i,ls)
c
c         termo:  +betah_p*(\mu_p,\lambda_p)
c
            eleffd(ncj,nci) = eleffd(ncj,nci) + betap*din*djn

c
         end do
c
 8111 continue
 1000 continue
 4000 continue
c
c   Condensacao
c
       call condtdg(elmbb,elmcb,elmdb,
     &            elfbb,eleffd,elresd,
     &            nepcon,nee)
c
c      computation of Dirichlet b.c. contribution
c
       call ztest(dl,nee,zerodl)
c
c       if(nel.eq.1) then
c           do i=1,nee
c               do j=1,nee
c                   write(45,*) i,j,eleffd(i,j),eleffd(j,i)
c               end do
c           end do
c       end if
c
      if(.not.zerodl)
     & call kdbc(eleffd,elresd,dl,nee)
c
c
c.... assemble element stifness matrix and force array into global
c        left-hand-side matrix and right-hand side vector
c
      call addlhs(alhs,eleffd,idiag,lm(1,1,nel),nee,diag)
c      call addnsl(alhs,dlhs,eleffd,idiag,lm(1,1,nel),nee,diag)
c
      call addrhs(brhs,elresd,lm(1,1,nel),nee)
c
  500 continue
      return
      end
c**** new **********************************************************************
      subroutine darcy3(ien   ,x     ,xl  ,
     &                 d     ,dl    ,mat   ,
     &                 matside,mcont,
     &                 det   ,shl   ,shg   ,
     &                 w     ,c     ,
     &                 grav  ,ipar  ,lado  ,
     &                 detc  ,shlc  ,shgc  ,
     &                 eleffd,elresd,shln  ,
     &                 shgn  ,wn    ,detn  ,
     &                 detb  ,shlb  ,shgb  ,
     &                 detpn ,shlpn ,shgpn ,
     &                 idside,xls   ,idlsd ,
     &                 dsfl  ,ddis  ,detp  ,
     &                 shlp  ,shgp  ,
     &                 elmbb ,
     &                 elmcb ,elmcc ,
     &                 elfbb ,elfcb ,
     &                 elfcc ,elmdb ,
c     
     &                 ideg  ,
c
     &                 shlpsd,shlcsd,
     &                 shgpsd,shgcsd,
c     
     &                 shlsd,
     &                 shgsd,
c
cc     &                 epres,
c     
     &                 numel ,neesq ,nen   ,
     &                 nsd   ,nesd  ,nint  ,
     &                 neg   ,nrowsh,ned   ,
     &                 nee   ,numnp ,ndof  ,
     &                 ncon  ,nencon,necon ,
     &                 neep  ,nints ,nnods ,
     &                 nenlad,npars ,nside ,
     &                 nenp  ,nodsp ,index ,
     &  	           nepcon,nface)
c
c
c.... program to calculate stifness matrix and force array for the
c        plane elasticity element and
c        assemble into the global left-hand-side matrix
c        and right-hand side vector
c
      implicit real*8 (a-h,o-z)
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      logical diag,quad,zerodl
      dimension elmbb(nepcon,*),elmcb(nepcon,*)
      dimension elmdb(nee,*),elmcc(nee,*)
      dimension elfbb(*),elfcb(*),elfcc(*)
      dimension ien(nen,*),x(nsd,*),xl(nesd,*),d(ndof,*),dl(ndof,*),
     &          mat(*),det(*),shl(nrowsh,nen,*),shg(nrowsh,nen,*),
     &          w(*),c(9,*),matside(*),mcont(*),
     &          grav(*),ipar(nodsp,*),lado(nside,*),detc(*),
     &          shlc(nrowsh,nencon,*),shgc(nrowsh,nencon,*),
     &          eleffd(nee,*),elresd(*)
      dimension shln(3,nnods,*),wn(*),detn(*),shgn(3,nnods,*),
     &          detb(*),shlb(3,nenlad,*),shgb(3,nenlad,*),
     &          detpn(*),shlpn(3,npars,*),shgpn(3,npars,*),
     &          idside(nside,*),xls(nesd,*),idlsd(*),
     &          dsfl(ncon,nencon,*),ddis(ned,nenp,*) 
      dimension dls(ndof,12),detp(*),shlp(3,nenp,*),shgp(3,nenp,*)
      dimension ue(6),duex(6),duey(6)
c
      dimension shlpsd(nrowsh,nenp,*),shlcsd(nrowsh,nencon,*),
     &          shgpsd(nrowsh,nenp,*),shgcsd(nrowsh,nencon,*)
c
      dimension ideg(2*ndof,*)
c
      dimension shlsd(nrowsh,nen,*),
     &          shgsd(nrowsh,nen,*)
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi

c
c      consistent matrix
c
      diag = .false.
      pi=4.d00*datan(1.d00)
      gf1=grav(1)
      gf2=grav(2)
      gf3=grav(3)      
c
c
      do 500 nel=1,numel
c
c      clear stiffness matrix and force array
c
      call clear(elmbb,nepcon*nepcon)
      call clear(elfbb,nepcon)
      call clear(elfcb,nee)
c
c      localize coordinates and Dirichlet b.c.
c
      call local(ien(1,nel),x,xl,nen,nsd,nesd)
      call local(ipar(1,nel),d,dl,nodsp,ndof,ndof)
c
      m = mat(nel)
      quad = .true.
      if (nen.eq.4.and.ien(3,nel).eq.ien(4,nel)) quad = .false.
c     
      call shgqs(xl,detc,shlc,shgc,nint,nel,neg,quad,nencon,shl,nen)
      call shgqs(xl,detp,shlp,shgp,nint,nel,neg,quad,nenp,shl,nen)
c
      nintb=nside*nints
!       call shgqsd(xl,shlpsd,shgpsd,nintb,nel,neg,nenp,shlsd,nen)
!       call shgqsd(xl,shlsd,shgsd,nintb,nel,neg,nen,shlsd,nen)
c
      call shgtqsd(ien(1,nel),xl,shlpsd,shgpsd,
     &     nside,nintb,nints,nel,neg,nenp,shlsd,nen)
      call shgtqsd(ien(1,nel),xl,shlcsd,shgcsd,
     &     nside,nintb,nints,nel,neg,nencon,shlsd,nen)

c
c....... form stiffness matrix
c
c      calculo de h - caracteristico para cada elemento
c
      if(nen.eq.3.) then
c
        h2=xl(1,2)*xl(2,3)+xl(1,1)*xl(2,2)+xl(1,3)*xl(2,1)
     &    -xl(1,1)*xl(2,3)-xl(1,2)*xl(2,1)-xl(1,3)*xl(2,2)
        h2=h2*pt5
        h=dsqrt(h2)
c
      else if (nen.eq.4) then
c
        h2=0
        do 100 l=1,2
          h2=h2+(xl(1,l)-xl(1,l+2))**2+(xl(2,l)-xl(2,l+2))**2
 100    continue 
        h=dsqrt(h2)/2.d00
        h2=h*h
c
      else 
        stop
      end if
c

c------------------------fonte Darcy tracador---------------------------
c
      if (nel.eq.1) then
        fd = 200.0d00
        
      else if (nel.eq.numel) then
        fd = -200.0d00
        
      else
        fd = 0.0d00
        
      end if
c-----------------------------------------------------------------------




c      set up material properties
c     
	  eps=c(3,m)  !
 	  delta1=c(4,m)
 	  delta2=c(5,m)
	  del1=c(6,m)
 	  del2=c(7,m)
      del3=c(8,m)
c
c     Matriz de condutividade anisotropica
c
c            | xka1  xka2 |
c K1 = alpha*|            | ,
c            | xka2  xka1 |
c
        alpha=c(9,m)
	  xka1=c(1,m)*alpha
	  xka2=c(2,m)*alpha
c
                xkl = c(1,m)*alpha
c
	    betap = xkl*delta1*h**delta2     !betah
c
c
c...  Matriz de Resistencia Hidraulica ({K}-1)
c
c               | a11  a12 |
c inv{K1} = A = |          | ,
c               | a21  a22 |
c
      deta=xka1*xka1-xka2*xka2
      a11 = xka1/deta
      a12 =-xka2/deta
      a22 = a11
      a21 = a12
c
c.....loop on integration points
c
      do 400 l=1,nint
      c1=detp(l)*w(l)
c      
      xx=0.d00
      yy=0.d00
c
      do 2323 i=1,nen
      xx=xx+shl(3,i,l)*xl(1,i)
      yy=yy+shl(3,i,l)*xl(2,i)
2323  continue
c
      pix=pi*xx
      piy=pi*yy
      sx=dsin(pix)
      sy=dsin(piy)
      cx=dcos(pix)
      cy=dcos(piy)
c
      pi2=pi*pi
c
c
c     fonte Darcy
c-----------------------------------------------------------------------
c      if (index.eq.1) then
c
c        if(m.eq.1) then
c            fd = 4.d00*pi2*sx*sy
c
c        else
c            fd = 2.d00*pi2*alpha*(2.d00*sx*sy - cx*cy)
c
c        end if
c
c      else if(index.eq.2) then
c
c      sx=dsin(xx)
c      sy=dsin(yy)
c      cx=dcos(xx)
c      cy=dcos(yy)
c      ex=dexp(xx)
c
c        if(m.eq.1) then
c            fd = alpha*xx*(2.d00*sy+cy) + sy
c
c        else
c            fd = -2.d00*alpha*ex*cy
c
c        end if
c
c      else if(index.eq.3) then
c            fd = cx*cy
c
c      end if
c-----------------------------------------------------------------------      
c
c      loop to compute volume integrals
c
      do 3000 j=1,nencon
      djx=shgc(1,j,l)*c1
      djy=shgc(2,j,l)*c1
      djn=shgc(3,j,l)*c1
c
c     curl(Av)
c
      rtj1=a11*djy - a21*djx
      rtj2=a12*djy - a22*djx
c
c.... source terms      
c
      naj=ncon*(j-1)
      naj1=naj+1
      naj2=naj+2
c
c     (f,div v)/k1
c
c   v1*f1
      elfbb(naj1)=elfbb(naj1) + del1*djx*fd/xkl
c   v2*f2
      elfbb(naj2)=elfbb(naj2) + del1*djy*fd/xkl
c
c
c.... element stiffness
c
      do 3003 i=1,nencon
      nai=ncon*(i-1)
      nai1=nai+1
      nai2=nai+2
c
      dix=shgc(1,i,l)
      diy=shgc(2,i,l)
      din=shgc(3,i,l)
c
c     curl(Au)
c
      rti1=a11*diy - a21*dix
      rti2=a12*diy - a22*dix
c
c
c termos:  (Av,u)
c
c
      elmbb(naj1,nai1) = elmbb(naj1,nai1)
     &                 + a11*djn*din
c
      elmbb(naj1,nai2) = elmbb(naj1,nai2)
     &                 + a12*djn*din
c
      elmbb(naj2,nai1) = elmbb(naj2,nai1)
     &                 + a21*djn*din
c
      elmbb(naj2,nai2) = elmbb(naj2,nai2) 
     &                 + a22*djn*din
c
c termos:  del3*K(Av,Au)
c
c
      elmbb(naj1,nai1) = elmbb(naj1,nai1)
     &                 + del3*djn*din*a11
c
      elmbb(naj1,nai2) = elmbb(naj1,nai2)
     &                 + del3*djn*din*a12
c
      elmbb(naj2,nai1) = elmbb(naj2,nai1)
     &                 + del3*djn*din*a21
c
      elmbb(naj2,nai2) = elmbb(naj2,nai2)
     &                 + del3*djn*din*a22
c
c
c termos:  del1*(div(v), div(u))/k1
c
c     v1*u1: 
      elmbb(naj1,nai1) = elmbb(naj1,nai1) 
     &                 + del1*djx*dix/xkl
c     v2*u1:
      elmbb(naj2,nai1) = elmbb(naj2,nai1)
     &                 + del1*djy*dix/xkl
c     v1*u2:
      elmbb(naj1,nai2) = elmbb(naj1,nai2)
     &                 + del1*djx*diy/xkl
c     v2*u2:
      elmbb(naj2,nai2) = elmbb(naj2,nai2) 
     &                 + del1*djy*diy/xkl
c
c	  
c
c termos: del2*(k1 rot(Av), rot(Au))
c
c     v1*u1: 
      elmbb(naj1,nai1) = elmbb(naj1,nai1) 
     &                 + del2*rtj1*rti1*xkl
c     v2*u1:
      elmbb(naj2,nai1) = elmbb(naj2,nai1)
     &                 + del2*rtj2*rti1*xkl
c     v1*u2:
      elmbb(naj1,nai2) = elmbb(naj1,nai2)
     &                 + del2*rtj1*rti2*xkl
c     v2*u2:
      elmbb(naj2,nai2) = elmbb(naj2,nai2) 
     &                 + del2*rtj2*rti2*xkl
c
 3003 continue
c
      do 3005 i=1,nenp
	  nbi1=i+necon
c       
      dix=shgp(1,i,l)
      diy=shgp(2,i,l)
      din=shgp(3,i,l)
c
c
c termos:  -(p, div(v)) 
c
c     p*v1: 
      elmbb(naj1,nbi1) = elmbb(naj1,nbi1) - din*djx 
c
c     p*v2: 
      elmbb(naj2,nbi1) = elmbb(naj2,nbi1) - din*djy
c
c
c termos:  del3*K*(gradp, Av)
c
c     p*v1: 
      elmbb(naj1,nbi1) = elmbb(naj1,nbi1) + del3*dix*djn
c
c     p*v2: 
      elmbb(naj2,nbi1) = elmbb(naj2,nbi1) + del3*diy*djn
c
c
 3005 continue
 3000 continue
c
c
c
      do 303 j=1,nenp
	  nbj1=j+necon
      djx=shgp(1,j,l)*c1
      djy=shgp(2,j,l)*c1
      djn=shgp(3,j,l)*c1
c
c     (f,q)
c
      elfbb(nbj1)=elfbb(nbj1) - fd*djn
c
      do 304 i=1,nencon
      nai=ncon*(i-1)
      nai1=nai+1
      nai2=nai+2
c
      dix=shgc(1,i,l)
      diy=shgc(2,i,l)
      din=shgc(3,i,l)
c
c
c termos:  -(div(u),q) 
c
c     u1*q: 
      elmbb(nbj1,nai1) = elmbb(nbj1,nai1) - djn*dix 
c
c
c     u2*q: 
      elmbb(nbj1,nai2) = elmbb(nbj1,nai2) - djn*diy
c
c
c termos:  del3*K*(Au,gradq)
c
c     u1*q: 
      elmbb(nbj1,nai1) = elmbb(nbj1,nai1) + del3*djx*din
c
c
c     u2*q: 
      elmbb(nbj1,nai2) = elmbb(nbj1,nai2) + del3*djy*din
c
c
  304 continue
      do 305 i=1,nenp
	  nbi1=i+necon     
      dix=shgp(1,i,l)
      diy=shgp(2,i,l)
      din=shgp(3,i,l)
c
c     termo: - eps*(p,q)
c
c     p*q: 
      elmbb(nbj1,nbi1) = elmbb(nbj1,nbi1) - eps*din*djn
c
c
c
c termos:  del3*K*(gradp,gradq)
c
      elmbb(nbj1,nbi1) = elmbb(nbj1,nbi1)
     &                 + del3*xka1*(dix*djx + diy*djy)
     &                 + del3*xka2*(diy*djx + dix*djy)
c
c
  305 continue
  303 continue
c
  400 continue
c
c......Symmetrization - boundary terms 
c
      do 4000 ns=1,nside
c
c-----localiza os parametros do lado ns
c
      do nn=1,npars
         nld = (ns-1)*npars + nn
         do jj=1,ndof
	  dls(jj,nn) = dl(jj,nld)
		 end do
      end do
c
c-----localiza os no's do lado ns
c
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,nel)
      nl2=ien(ns2,nel)
c      
      if(nl2.gt.nl1) then 
	  sign = 1.d00
        do  nn=1,nenlad
           idlsd(nn)=idside(ns,nn)
	  end do
	else
	  sign = -1.d00
	  idlsd(1) = idside(ns,2)
	  idlsd(2) = idside(ns,1)
	  id3 = nenlad-2
	  if(id3.gt.0) then
	    do il=3,nenlad
	      idlsd(il) = idside(ns,nenlad+3-il)
	    end do 
	  end if
      end if 
c
c    localize the coordinates of the side nodes 
c
      do 2000 nn=1,nenlad
      nl=idlsd(nn)
      xls(1,nn)=xl(1,nl)
      xls(2,nn)=xl(2,nl)
 2000 continue
c
      call oneshgp(xls,detn,shlb,shln,shgn,
     &             nenlad,nnods,nints,nesd,ns,nel,neg)      
      call oneshgp(xls,detpn,shlb,shlpn,shgpn,
     &             nenlad,npars,nints,nesd,ns,nel,neg)
c
c.....compute boundary integral - symmetrization
c
      do 1000 ls = 1,nints
      lb = (ns-1)*nints + ls
c
c    valores dos parametros do multiplicador two degree of freedom
c
        dhs = 0.d00
          do i=1,npars
	      dhs = dhs + dls(1,i)*shgpn(2,i,ls)
          end do
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c
c
c TESTE: valores exatos do multiplicador
c

         call uexafx(x1,x2,ue,duex,duey,alpha,m,index)
c
c
       dhse=ue(3)            ! p
c
c       dhs = dhse
c
       do 1110 j=1,nencon
          naj = ncon*(j-1)
          naj1=naj+1
          naj2=naj+2
          djn=shgcsd(3,j,lb)*detn(ls)*wn(ls)
          djx=shgcsd(1,j,lb)*detn(ls)*wn(ls)
          djy=shgcsd(2,j,lb)*detn(ls)*wn(ls)
          gjn=djx*xn1 + djy*xn2
c
c         termo: -(p_exato,v.n)
c
          elfbb(naj1)= elfbb(naj1) - dhs*djn*xn1
c
          elfbb(naj2)= elfbb(naj2) - dhs*djn*xn2
c
c
 1110 continue
c
c    jump terms
c
       do 9100 j=1,nenp
          nbj1=j+necon
          djn=shgpsd(3,j,lb)*detn(ls)*wn(ls)
          djx=shgpsd(1,j,lb)*detn(ls)*wn(ls)
          djy=shgpsd(2,j,lb)*detn(ls)*wn(ls)
          gjn=djx*xn1 + djy*xn2
          djn1=djn*xn1
          djn2=djn*xn2
c
c         termo:  betap*(q,p_exato)
c
          elfbb(nbj1) = elfbb(nbj1) + betap*dhs*djn
c
       do 9200 i=1,nenp
          nbi1=i+necon
          din=shgpsd(3,i,lb)
          dix=shgpsd(1,i,lb)
          diy=shgpsd(2,i,lb)
          gin=dix*xn1 + diy*xn2
          din1=din*xn1
          din2=din*xn2
c
c         termo:  betap*(p,q)
c
          elmbb(nbj1,nbi1) = elmbb(nbj1,nbi1) + betap*din*djn
c 
 9200 continue
 9100 continue
c
 1000 continue
 4000 continue
c
c	Local TLDG-solution (potencial)  
c
      call solvetdg(elmbb,elfbb,nepcon)
c          
c      valores nodais descontinuos
c
        jk=0
        do 9171 j=1,nencon
        do 9272 jj=1,ncon
        jk=jk+1
        dsfl(jj,j,nel)=elfbb(jk)
 9272   continue
 9171   continue
c
      do j=1,nenp
	  jn=j+necon
	  ddis(1,j,nel) = elfbb(jn)
      end do
c
c
  500 continue
c


c      write(*,*)'gera saida'
c-----------------------------------------------------------------------    
      do i = 1, numel
        do j = 1, nen
          do k =1, nsd
            write(*,*)dsfl(k, j, i)
c            write(*,*)i, j, k  
          end do
        end do
      end do
c-----------------------------------------------------------------------

      return
      end
      
c**** new **********************************************************************
      subroutine local(ien,x,xl,nen,nrowx,nrowxl)
c
c.... program to localize a global array
c
c        note: it is assumed nrowxl.le.nrowx
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension ien(*),x(nrowx,*),xl(nrowxl,*)
c
      do 200 j=1,nen
      node = ien(j)
c
      do 100 i=1,nrowxl
      xl(i,j)= x(i,node)
  100 continue
c
  200 continue
c
      return
      end
    
c***********************************************************************
c             miscelânea
c***********************************************************************      

c**** new **********************************************************************
      block data
c
c.... program to define output labels and numerical constants
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /labels/ labeld(3),label1(16),label2(3)
c
c        labeld(3)  = displacement, velocity and acceleration labels
c        label1(16) = output labels for element-type 1
c        label2(3)  = output labels for element-type 2
c
c.... note: add label arrays for any additional elements        
c
      data   zero,pt1667,pt25,pt5
     &      /0.0d0,0.1666666666666667d0,0.25d0,0.5d0/,
     &       one,two,three,four,five,six
     &      /1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0/
c
      end
c**** new ***************************************************************
      subroutine echo 
c 
c.... program to echo input data 
c 
      character*4 ia(20) 
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c 
      read(iin,1000) iech 
      if (iech.eq.0) return 
c 
      write(iecho,2000) iech 
      backspace iin 
c 
      do 100 i=1,100000 
      read(iin,3000,end=200) ia 
      if (mod(i,50).eq.1) write(iecho,4000) 
      write(iecho,5000) ia 
  100 continue 
c 
  200 continue 
      rewind iin 
      read(iin,1000) iech 
c 
      return 
c 
 1000 format(16i10) 
 2000 format(///,' i n p u t   d a t a   f i l e               ',  //5x, 
     &' echo print code . . . . . . . . . . . . . . (iecho ) = ',i10//5x, 
     &'    eq. 0, no echo of input data                        ',   /5x, 
     &'    eq. 1, echo input data                              ',   ///) 
 3000 format(20a4) 
 4000 format(' ',8('123456789*'),//) 
 5000 format(' ',20a4) 
      end 
c************************************************************************** 
      function mpoint(name,ndim1,ndim2,ndim3,ipr) 
c 
c.... program to calculate storage pointer 
c 
      character*4 name(2) 
      common /bpoint/ mfirst,mlast,ilast,mtot,iprec
c 
      mpoint = mfirst 
      if ( iprec.eq.2 .and. mod(mpoint,2).eq.0 ) mpoint = mpoint + 1 
      call dctnry(name,ndim1,ndim2,ndim3,mpoint,ipr,mlast,ilast) 
      mfirst = mpoint + ndim1*max0(1,ndim2)*max0(1,ndim3)*ipr 
      if (mfirst.ge.mlast) call serror(name,mfirst-mlast) 
c 
      return 
      end 
c************************************************************************** 
      subroutine dctnry(name,ndim1,ndim2,ndim3,mpoint,ipr,mlast,ilast) 
c 
c.... program to store pointer information in dictionary 
c 
      character*4 name(2) 
      character*4 ia 
      common na(1) 
      common /dictn/ ia(1) 
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c 
      mlast = mlast - 5 
      ia(ilast+1) = name(1) 
      ia(ilast+2) = name(2) 
      na(mlast+1) = mpoint 
      na(mlast+2) = ndim1 
      na(mlast+3) = ndim2 
      na(mlast+4) = ndim3 
      na(mlast+5) = ipr 
      ilast = ilast + 2 
c 
      return 
      end 
c************************************************************************
      subroutine serror(name,i) 
c 
c.... program to print error message if available storage is exceeded 
c 
      character*4 name(2) 
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c 
      call prtdc 
      write(*,1000) i,name
      pause 
      stop 
c 
 1000 format(1x,5('*'),'storage exceeded by ',i10, 
     &' words in attempting to store array ',2a4) 
      end 
c**** new ********************************************************************** 
      subroutine prtdc 
c 
c.... program to print memory-pointer dictionary 
c 
      common /bpoint/ mfirst,mlast,ilast,mtot,iprec 
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
      character *4 ia 
      common na(1) 
      common /dictn/ ia(1) 
c 
      n = (mtot-mlast)/5 
      j = mtot + 1 
c 
      k = 1 
      do 100 i=1,n 
      if (mod(i,50).eq.1) write(iout,1000) 
      j = j - 5 
      call prtdc1(i,ia(k),na(j),na(j+1),na(j+2),na(j+3),na(j+4)) 
      k = k + 2 
 100  continue 
c 
      return 
c 
 1000 format(///, 
     &' d y n a m i c   s t o r a g e    a l l o c a t i o n', 
     &'   i n f o r m a t i o n '// 
     &  12x,'array no.',5x,'array',8x,'address',6x,'dim1',6x,'dim2', 
     &  6x, 'dim3',6x,'prec.'/) 
c 
      end 
c**** new ********************************************************************** 
      subroutine prtdc1(i,iname,iadd,ndim1,ndim2,ndim3,ipr) 
c 
c.... program to print memory-pointer information for an array 
c 
      character *4 iname(2) 
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c 
      if (i.eq.1) neg = 1 
      if (iname(1).eq.'npar') then 
        write (iout,1000) neg 
        neg = neg + 1 
      endif 
      write(iout,2000) i,iname,iadd,ndim1,ndim2,ndim3,ipr 
c 
      return 
c 
 1000 format(/14x,'*****',7x,'begin element group number',i10/' ') 
 2000 format(14x,i10,7x,2a4,1x,6i10) 
      end 
      
c**** new *********************************************************************
      subroutine imove(ia,ib,n)
c
c.... program to move an integer array 
c
      dimension ia(*),ib(*)
c
      do 100 i=1,n
      ia(i)=ib(i)
  100 continue

      return
      end
      
c**** new **********************************************************************
      subroutine move(a,b,n)
c
c.... program to move a floating-point array
c
      implicit real*8(a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension a(*),b(*)
c
      do 100 i=1,n
      a(i) = b(i)
  100 continue
c
      return
      end
c**** new **********************************************************************
      subroutine iclear(ia,m)
c
c.... program to clear an integer array
c
      dimension ia(*)
c
      do 100 i=1,m
      ia(i) = 0
  100 continue
c
      return
      end
c**** new **********************************************************************
      subroutine clear(a,m)
c
c.... program to clear a floating-point array
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension a(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      do 100 i=1,m
      a(i) = zero
  100 continue
c
      return
      end
c**** new **********************************************************************
      subroutine ztest(a,n,lzero)
c
c.... program to determine if an array contains only zero entries
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension a(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      logical lzero
c
      lzero = .true.
c
      do 100 i=1,n
      if (a(i).ne.zero) then
         lzero = .false.
         return
      endif
  100 continue
c
      end     
c
c**** new **********************************************************************
      subroutine multab(a,b,c,ma,mb,mc,l,m,n,iopt)
c
c.... program to multiply two matrices
c
c        l = range of dot-product index
c        m = number of active rows in c
c        n = number of active columns in c
c
      implicit real*8 (a-h,o-z)
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      dimension a(ma,*),b(mb,*),c(mc,*)
c
      go to (1000,2000,3000,4000),iopt
c
c.... iopt = 1, c(i,j) = a(i,k)*b(k,j) , (c = a * b)
c
 1000 do 1200 i=1,m
c
      do 1100 j=1,n
      c(i,j) = rcdot(a(i,1),b(1,j),ma,l)
 1100 continue
c
 1200 continue
      return
c                                            t
c.... iopt = 2, c(i,j) = a(k,i)*b(k,j) (c = a  * b)
c
 2000 do 2200 i=1,m
c
      do 2100 j=1,n
      c(i,j) = coldot(a(1,i),b(1,j),l)
 2100 continue
c
 2200 continue
      return
c                                                t
c.... iopt = 3, c(i,j) = a(i,k)*b(j,k) (c = a * b )
c
 3000 do 3200 i=1,m
c
      do 3100 j=1,n
      c(i,j) = rowdot(a(i,1),b(j,1),ma,mb,l)
 3100 continue
c
 3200 continue
      return
c                                            t    t
c.... iopt = 4, c(i,j) = a(k,i)*b(j,k) (c = a  * b )
c
 4000 do 4200 i=1,m
c
      do 4100 j=1,n
      c(i,j) = rcdot(b(j,1),a(1,i),mb,l)
 4100 continue
c
 4200 continue
c
      return
      end
c**** new **********************************************************************
      function rowdot(a,b,ma,mb,n)
c
c.... program to compute the dot product of vectors stored row-wise
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension a(ma,*),b(mb,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      rowdot = zero
c
      do 100 i=1,n
      rowdot = rowdot + a(1,i)*b(1,i)
!       write(58,*)a(1,i),'*',b(1,i),'=',rowdot
  100 continue
c
      return
      end    
c**** new **********************************************************************
      function rcdot(a,b,ma,n)
c
c.... program to compute the dot product of a vector stored row-wise
c        with a vector stored column-wise
c
      implicit real*8 (a-h,o-z)
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      dimension a(ma,*),b(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      rcdot = zero
c
      do 100 i=1,n
      rcdot = rcdot + a(1,i)*b(i)
  100 continue
c
      return
      end
c**** new **********************************************************************
      function coldot(a,b,n)
c
c.... program to compute the dot product of vectors stored column-wise
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension a(*),b(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      coldot = zero
c
      do 100 i=1,n
      coldot = coldot + a(i)*b(i)
  100 continue
c
      return
      end

c***********************************************************************
c             geometria
c***********************************************************************

      subroutine coord(x,nsd,numnp,iprtin)
c
c.... program to read, generate and write coordinate data
c
c        x(nsd,numnp) = coordinate array
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension x(nsd,*)
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
      call genfl(x,nsd)
c
      if (iprtin.eq.1) return
c
      do 100 n=1,numnp
      if (mod(n,50).eq.1) write(iecho,1000) (i,i=1,nsd)
      write(iecho,2000) n,(x(i,n),i=1,nsd)   
  100 continue
c
      return
c
 1000 format(///,' n o d a l   c o o r d i n a t e   d a t a '///5x,
     &' node no.',3(13x,' x',i1,' ',:)//)
 2000 format(6x,i10,10x,3(1pe15.8,2x))
      end
c******************************************************************************
      subroutine genfl(a,nra)                                           
c                                                                       
c.... program to read and generate floating-point nodal data            
c                                                                       
c         a       = input array                                         
c         nra     = number of rows in a (le.6)                          
c         n       = node number                                         
c         numgp   = number of generation points                         
c         ninc(i) = number of increments for direction i                
c         inc(i)  = increment for direction i                           
c                                                                       
      implicit real*8(a-h,o-z)                                          
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      dimension a(nra,*)                                                
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
      common /genflc/ temp(6,20),n,numgp,ninc(3),inc(3)                 
c                                                                       
  100 continue                                                          
      read(iin,1000) n,numgp,(temp(i,1),i=1,nra)                        
      if (n.eq.0) return                                                
      call move(a(1,n),temp,nra)                                       
      if (numgp.ne.0) then
         do 200 j=2,numgp                                                  
c                                                                       
         read(iin,1000) m,mgen,(temp(i,j),i=1,nra)                         
         if (mgen.ne.0) call move(temp(1,j),a(1,m),nra) 
c                                                                       
  200    continue                                                          
         read(iin,2000) (ninc(i),inc(i),i=1,3)                             
         call genfl1(a,nra)                                                
      endif
      go to 100                                                         
c                                                                       
 1000 format(2i10,6f10.0)                                                
 2000 format(16i10)                                                      
c                                                                       
      end                                                               
c******************************************************************************
      subroutine genfl1(a,nra)                                          
c                                                                       
c.... program to generate floating-point nodal data 
c        via isoparametric interpolation         
c                                                                       
c         iopt = 1, generation along a line                             
c              = 2, generation over a surface                           
c              = 3, generation within a volume                            
c                                                                       
      implicit real*8(a-h,o-z)                                          
c                                                                       
c.... remove above card for single-precision operation                  
c                                                                       
      dimension a(nra,*),sh(20)                                         
      common /genflc/ temp(6,20),n,numgp,ninc(3),inc(3)                 
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      iopt = 3                                                          
      if (ninc(3).eq.0) iopt = 2                                        
      if (ninc(2).eq.0) iopt = 1                                        
c                                                                       
      dr = zero                                         
      ds = zero                                                        
      dt = zero                                   
c                                                                       
      if (ninc(1).ne.0) dr = two/ninc(1)                                
      if (ninc(2).ne.0) ds = two/ninc(2)                                
      if (ninc(3).ne.0) dt = two/ninc(3)                                
c                                                                       
      ii = ninc(1)+1                                                    
      jj = ninc(2)+1                                                    
      kk = ninc(3)+1                                                    
c                                                                       
      ni = n                                                            
      nj = n                                                            
      nk = n                                                            
c                                                                       
      t = -one                                                          
      do 300 k=1,kk                                                     
c
      s = -one                                                          
      do 200 j=1,jj                                                     
c
      r = -one                                                          
      do 100 i=1,ii                                                     
c                                                                       
      call gensh(r,s,t,sh,numgp,iopt)                                   
      call multab(temp,sh,a(1,ni),6,20,nra,numgp,nra,1,1)               
      ni = ni + inc(1)                                                      
      r = r + dr                                                            
  100 continue                                                          
c                                                                       
      nj = nj + inc(2)                                                      
      ni = nj                                                             
      s = s + ds                                                            
  200 continue                                                          
c                                                                       
      nk = nk + inc(3)                                                      
      ni = nk                                                             
      t = t + dt                                                            
  300 continue                                                          
c                                                                       
      return                                                            
      end                                                               
c****************************************************************************** 
    
      subroutine gensh(r,s,t,sh,numgp,iopt)                             
c                                                                       
c.... program to call shape function routines         
c        for isoparametric generation         
c                                                                       
      implicit real*8(a-h,o-z)                                      
c                                                                       
c.... modify above card for single-precision operation               
c                                                                       
      dimension sh(*)                                                   
c                                                                       
      go to (100,200,300),iopt                                                
c                                                                       
  100 call gensh1(r,sh,numgp)                                           
      return                                                            
c                                                                       
  200 call gensh2(r,s,sh,numgp)                                         
      return                                                            
c                                                                       
  300 call gensh3(r,s,t,sh,numgp)                                       
      return                                                            
c                                                                       
      end                                                               
c******************************************************************************
      subroutine gensh1(r,sh,n)                                         
c                                                                       
c.... program to compute 1d shape functions           
c        for isoparametric generation                     
c                                                                       
      implicit real*8(a-h,o-z)                                          
c                                                                       
c.... modify above card(s) for single-precision operation               
c                                                                       
      dimension sh(*)                                                   
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c                                                                       
      sh(2) = pt5*r                                                       
      sh(1) = pt5 - sh(2)                                                   
      sh(2) = pt5 + sh(2)                                                   
      if (n.eq.3) then
         sh(3) = one - r*r                                                     
         sh(1) = sh(1) - pt5*sh(3)                                             
         sh(2) = sh(2) - pt5*sh(3)                                             
      endif
c                                                                       
      return                                                            
      end                                                               
c******************************************************************************
      subroutine gensh2(r,s,sh,n)                                       
c                                                                       
c.... program to compute 2d shape functions 
c        for isoparametric generation    
c                                                                       
      implicit real*8(a-h,o-z)                                          
c                                                                       
c.... modify above card for single-precision operation               
c                                                                       
      dimension sh(*)                                      
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      r2 = pt5*r                                                          
      r1 = pt5 - r2                                                         
      r2 = pt5 + r2                                                         
      s2 = pt5*s                                                          
      s1 = pt5 - s2                                                         
      s2 = pt5 + s2                                                         
      sh(1) = r1*s1                                                       
      sh(2) = r2*s1                                                       
      sh(3) = r2*s2                                                       
      sh(4) = r1*s2                                                       
      if (n.eq.4) return                                                
c                                                                       
      r3 = one - r*r                                                        
      s3 = one - s*s                                                        
      sh(5) = r3*s1                                                       
      sh(6) = s3*r2                                                       
      sh(7) = r3*s2                                                       
      sh(8) = s3*r1                                                       
      sh(1) = sh(1) - pt5*(sh(5) + sh(8))
      sh(2) = sh(2) - pt5*(sh(6) + sh(5))
      sh(3) = sh(3) - pt5*(sh(7) + sh(6))
      sh(4) = sh(4) - pt5*(sh(8) + sh(7))
c                                                                       
      return                                                            
      end                                                               
c******************************************************************************
      subroutine gensh3(r,s,t,sh,n)                                     
c                                                                       
c.... program to compute 3d shape functions            
c        for isoparametric generation   
c                                                                       
      implicit real*8(a-h,o-z)                                          
c                                                                       
c.... modify above card for single-precision operation               
c                                                                       
      dimension sh(*)                                                   
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c                                                                       
      r2 = pt5*r
      r1 = pt5 - r2                                                         
      r2 = pt5 + r2                                                         
      s2 = pt5*s
      s1 = pt5 - s2                                                         
      s2 = pt5 + s2                                                         
      t2 = pt5*t
      t1 = pt5 - t2                                                         
      t2 = pt5 + t2                                                         
c                                                                       
      rs1 = r1*s1                                                         
      rs2 = r2*s1                                                         
      rs3 = r2*s2                                                         
      rs4 = r1*s2                                                         
      sh(1) = rs1*t1                                                      
      sh(2) = rs2*t1                                                      
      sh(3) = rs3*t1                                                      
      sh(4) = rs4*t1                                                      
      sh(5) = rs1*t2                                                      
      sh(6) = rs2*t2                                                      
      sh(7) = rs3*t2                                                      
      sh(8) = rs4*t2                                                      
      if (n.eq.8) return                                                 
c                                                                       
      r3 = one - r*r                                                        
      s3 = one - s*s                                                        
      t3 = one - t*t                                                        
      sh(17) = t3*rs1                                                     
      sh(18) = t3*rs2                                                     
      sh(19) = t3*rs3                                                     
      sh(20) = t3*rs4                                                     
      rs1 = r3*s1                                                         
      rs2 = s3*r2                                                         
      rs3 = r3*s2                                                         
      rs4 = s3*r1                                                         
      sh( 9) = rs1*t1                                                     
      sh(10) = rs2*t1                                                     
      sh(11) = rs3*t1                                                     
      sh(12) = rs4*t1                                                     
      sh(13) = rs1*t2                                                     
      sh(14) = rs2*t2                                                     
      sh(15) = rs3*t2                                                     
      sh(16) = rs4*t2                                                     
c                                                                       
      sh(1) = sh(1) - pt5*(sh( 9) + sh(12) + sh(17))
      sh(2) = sh(2) - pt5*(sh( 9) + sh(10) + sh(18))
      sh(3) = sh(3) - pt5*(sh(10) + sh(11) + sh(19))
      sh(4) = sh(4) - pt5*(sh(11) + sh(12) + sh(20))
      sh(5) = sh(5) - pt5*(sh(13) + sh(16) + sh(17))
      sh(6) = sh(6) - pt5*(sh(13) + sh(14) + sh(18))
      sh(7) = sh(7) - pt5*(sh(14) + sh(15) + sh(19))
      sh(8) = sh(8) - pt5*(sh(15) + sh(16) + sh(20))
c                                                                       
      return                                                            
      end                                                               
c***********************************************************************
c     boundary
c***********************************************************************
c**** new **********************************************************************
      subroutine bc(id,ndof,numnp,neq,iprtin)
c
c.... program to read, generate and write boundary condition data
c        and establish equation numbers
c
      dimension id(ndof,*)
c
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
      logical pflag
c
      call iclear(id,ndof*numnp)
      call igen(id,ndof)
c
      if (iprtin.eq.0) then
         nn=0
         do 200 n=1,numnp
         pflag = .false.
c
         do 100 i=1,ndof
         if (id(i,n).ne.0) pflag = .true.
  100    continue
c
         if (pflag) then      
            nn = nn + 1
            if (mod(nn,50).eq.1) write(iecho,1000) (i,i=1,ndof)
            write(iecho,2000) n,(id(i,n),i=1,ndof)
         endif
  200    continue
      endif
c
c.... establish equation numbers
c
      neq = 0
c
      do 400 n=1,numnp
c
      do 300 i=1,ndof
      if (id(i,n).eq.0) then
         neq = neq + 1
         id(i,n) = neq
      else
         id(i,n) = 1 - id(i,n)
      endif
c
  300 continue
c
  400 continue
c
      return
c
 1000 format(///,' n o d a l   b o u n d a r y   c o n d i t i o n  c o
     & d e s'///
     & 5x,' node no.',3x,6(6x,'dof',i1:)//)
 2000 format(6x,i10,5x,6(5x,i10))
c
      end
c***********************************************************************
      subroutine bcedge(ideg,id,npar,nedge,ndof,numnp,neq,iprtin)
c
c.... program to read, generate and write boundary condition data
c        and establish equation numbers
c
      dimension id(ndof,*),ideg(2*ndof,*)
c
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
      logical pflag
c
      call iclear(ideg,2*ndof*nedge)
      call iclear(id,ndof*numnp)
      call igen(ideg,2*ndof)
c
      if (iprtin.eq.0) then
         nn=0
         do 200 n=1,nedge
         pflag = .false.
c
         do 100 i=1,2*ndof
         if (ideg(i,n).ne.0) pflag = .true.
  100    continue
c
         if (pflag) then      
            nn = nn + 1
            if (mod(nn,50).eq.1) write(iecho,1000) (i,i=1,ndof)
            write(iecho,2000) n,(ideg(i,n),i=1,2*ndof)
         endif
  200    continue
      endif
c
c    id - prescribed dof
c
        kk=0
        do n=1,nedge
	     do is=1,npar
	       kk=kk+1
             do j=1,ndof
	         id(j,kk) = ideg(j,n)
	       end do
	     end do
	  end do
c
      if (iprtin.eq.0) then
         nn=0
         do 220 n=1,numnp
         pflag = .false.
c
         do 110 i=1,ndof
         if (id(i,n).ne.0) pflag = .true.
  110    continue
c
         if (pflag) then      
            nn = nn + 1
            if (mod(nn,50).eq.1) write(iecho,1100) (i,i=1,ndof)
            write(iecho,2000) n,(id(i,n),i=1,ndof)
         endif
  220    continue
      endif
c
c.... establish equation numbers
c
      neq = 0
c
      write(iecho,*)'Matriz ID:'
      do 400 n=1,numnp
c
      do 300 i=1,ndof
      if (id(i,n).eq.0) then
         neq = neq + 1
         id(i,n) = neq
      else
         id(i,n) = 1 - id(i,n)
      endif
c
  300 continue
c
      write(iecho,2000) n,(id(i,n),i=1,ndof)
  400 continue
c
      write(iecho,*)'neq = ',neq
      return
c
 1000 format(//, 10x, ' e d g e   b o u n d a r y   c o n d i t i o n  c o
     & d e s'///
     & 5x,'   node no.',3x,6(13x,'dof',i1:)//)
 1100 format(//, 10x,' n o d a l   b o u n d a r y   c o n d i t i o n  c o
     & d e s'///
     & 5x,'   node no.',3x,6(13x,'dof',i1:)//)
 2000 format(6x,i10,5x,6(5x,i10))
c
      end
c**** new **********************************************************************
      subroutine kdbc(eleffm,elresf,dl,nee)
c
c.... program to adjust load vector for prescribed displacement
c     boundary condition
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension eleffm(nee,*),elresf(*),dl(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      do 200 j=1,nee
c
      val=dl(j)
      if(val.eq.zero) go to 200
c
      do 100 i=1,nee
      elresf(i)=elresf(i)-eleffm(i,j)*val
100   continue
c
200   continue
c
      return
      end
c**** new **********************************************************************
      subroutine drchbc(shgpn,shlb,detpn,wn,
     &            xls,dls,pi,eps,alpha,ndgs,
     &            nints,nenlad,npars,ned,igl,m,index)
c                                                                        
c     Dirichlet boundary conditions                                             
c                                                                        
      implicit real*8(a-h,o-z)                                  
      dimension shlb(3,nenlad,*),shgpn(3,npars,*),xls(2,*)
      dimension detpn(*),wn(*),dls(ned,*)
      dimension aa(10,10),bb(10)
      dimension ue(6),duex(6),duey(6)
c

      do i=1,npars
	 bb(i) = 0.d00
	do j=1,npars
	 aa(i,j) = 0.d00
      end do
	end do
c
      do 1000 ls=1,nints
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
        end do
c
c    valor exato do multiplicador
c

         call uexafx(x1,x2,ue,duex,duey,alpha,m,index)
c
c
       dhs=ue(3)            ! p
c
c
        do 1100 j=1,npars
c
         djn=shgpn(2,j,ls)*detpn(ls)*wn(ls)
c
c    source term
c
         bb(j) = bb(j) + dhs*djn
c
        do i=1,npars
c
          din=shgpn(2,i,ls)
c
c     L2-projection at the boundary (side)
c             
          aa(i,j) = aa(i,j) + din*djn
c
        end do
c
 1100 continue
 1000 continue
c
      if(npars.gt.1) call invmb(aa,10,npars)
      if(npars.eq.1) aa(1,1) = 1.d00/aa(1,1)
c
	do i=1,npars
	  dls(igl,i) = 0.d00
	do j=1,npars      
	  dls(igl,i) = dls(igl,i) + aa(i,j)*bb(j)
	end do
	end do
      return                                                             
      end                                                                
c**** new **********************************************************************
c**** new **********************************************************************
      subroutine neuhbc(shgpn,shlb,detpn,wn,idlsd,
     &            xls,elresd,pi,eps,sign,ned,alpha,m,
     &            nints,nenlad,ns,npars,index)
c                                                                        
c     Neumann boundary conditions                                             
c                                                                        
      implicit real*8(a-h,o-z)                                  
      dimension shlb(3,nenlad,*),shgpn(3,npars,*),xls(2,*)
	  dimension detpn(*),wn(*),elresd(*),idlsd(*)
      dimension dhs(3),ue(6),duex(6),duey(6)
c
      do 1000 ls=1,nints
c
c     geometria
c
        un =0.d00
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c
         call uexafx(x1,x2,ue,duex,duey,alpha,m,index)
c
       dhs(1)=ue(1)            ! u1
       dhs(2)=ue(2)            ! u2
       dhs(3)=ue(3)            ! p
c

c	   un = dhs(1)*xn1+dhs(2)*xn2

c-----------------------------------------------------------------------
       un = 0.d00
c-----------------------------------------------------------------------
c
        do 1100 j=1,npars
          ncj = (ns-1)*npars  + j
c
          djn=shgpn(2,j,ls)*detpn(ls)*wn(ls)

c
c    source term
c        
          elresd(ncj) = elresd(ncj) + un*djn
c
 1100 continue
 1000 continue
c
      return                                                             
      end                                                                
c**** new **********************************************************************
c***********************************************************************
c       conectivity
c***********************************************************************

      subroutine igen(ia,m)
c
c.... program to read and generate integer nodal data
c
c        ia = input array
c         m = number of rows in ia
c         n = node number
c        ne = end node in generation sequence
c        ng = generation increment
c
      dimension ia(m,*),ib(13)
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
  100 continue
      read(iin,1000) n,ne,ng,(ib(i),i=1,m)
      if (n.eq.0) return
      if (ng.eq.0) then
         ne = n
         ng = 1
      else
         ne = ne - mod(ne-n,ng)
      endif
c
      do 200 i=n,ne,ng
      call imove(ia(1,i),ib,m)
  200 continue
c
      go to 100
c
 1000 format(16i10)
      end
c******************************************************************************
      subroutine genel(ien,mat,nen)                                     
c                                                                       
c.... program to read and generate element node and material numbers    
c                                                                       
c         ien(nen,numel) = element node numbers                         
c         mat(numel)     = element material numbers                     
c         nen            = number of element nodes (le.27)              
c         n              = element number                               
c         ng             = generation parameter                         
c         nel(i)         = number of elements in direction i            
c         incel(i)       = element number increment for direction i     
c         inc(i)         = node number increment for direction i        
c                                                                       
      dimension ien(nen,*),mat(*),itemp(27)                             
      common /genelc/ n,nel(3),incel(3),inc(3)                          
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
  100 continue
      read(iin,1000) n,m,(itemp(i),i=1,nen),ng
c      write(*,*) n,m,(itemp(i),i=1,nen),ng
      if (n.eq.0) return
      call imove(ien(1,n),itemp,nen)                                    
      mat(n)=m                                                          
      if (ng.ne.0) then
c                                                                       
c....... generate data                                                     
c                                                                       
       read(iin,1000) (nel(i),incel(i),inc(i),i=1,3)                     
c       write(*,*) (nel(i),incel(i),inc(i),i=1,3)
       call genel1(ien,mat,nen)
      endif
      go to 100                                                         
c                                                                       
 1000 format(20i10)                                             
c                                                                       
      end                                                               
c******************************************************************************
      subroutine genel1(ien,mat,nen)                                    
c                                                                       
c.... program to generate element node and material numbers             
c                                                                       
      dimension ien(nen,*),mat(*)                                       
      common /genelc/ n,nel(3),incel(3),inc(3)                          
c
 1000 format(16i10)                                             
c                                                                       
c.... set defaults                                                      
c                                                                       
      call geneld                                                       
c                                                                       
c.... generation algorithm                                              
c
      ie = n                                                            
      je = n                                                            
      ke = n                                                            
c                                                                      
      ii = nel(1)                                                       
      jj = nel(2)                                                       
      kk = nel(3)                                                       
c
c                                                                       
      do 300 k=1,kk                                                     
c
      do 200 j=1,jj                                                     
c
      do 100 i=1,ii                                                     
c
c                                                                       
      if (i.ne.ii) then
         le = ie            	                                                 
         ie = le + incel(1)                                                
         call geneli(ien(1,ie),ien(1,le),inc(1),nen)  
         mat(ie) = mat(le)                                                 
      endif
  100 continue                                                          
c                                                                       
      if (j.ne.jj) then
         le = je                                                           
         je = le + incel(2)                                                
         call geneli(ien(1,je),ien(1,le),inc(2),nen)                       
         mat(je) = mat(le)                                                 
         ie = je                                                           
      endif
  200 continue                                                          
c                                                                             
      if (k.ne.kk) then
         le = ke                                                           
         ke = le + incel(3)                                                
         call geneli(ien(1,ke),ien(1,le),inc(3),nen)                       
         mat(ke) = mat(le)                                                 
         ie = ke                                                           
      endif
  300 continue                                                          
c                                                                       
      return                                                            
      end                                                               
c******************************************************************************
      subroutine geneld                                                 
c                                                                       
c.... program to set defaults for element node       
c        and material number generation                              
c                                                                       
      common /genelc/ n,nel(3),incel(3),inc(3)                          
c                                                                       
      if (nel(1).eq.0) nel(1) = 1                                       
      if (nel(2).eq.0) nel(2) = 1                                       
      if (nel(3).eq.0) nel(3) = 1                                       
c                                                                       
      if (incel(1).eq.0) incel(1) = 1                                   
      if (incel(2).eq.0) incel(2) = nel(1)                              
      if (incel(3).eq.0) incel(3) = nel(1)*nel(2)                       
c                                                                       
      if (inc(1).eq.0) inc(1) = 1                                       
      if (inc(2).eq.0) inc(2) = (1+nel(1))*inc(1)                       
      if (inc(3).eq.0) inc(3) = (1+nel(2))*inc(2)                       
c                                                                       
      return                                                            
      end                                                               
c******************************************************************************
      subroutine geneli(ien2,ien1,inc,nen)                              
c                                                                       
c.... program to increment element node numbers                         
c                                                                       
      dimension ien1(*),ien2(*)                                         
c                                                                       
      do 100 i=1,nen                                                    
      if (ien1(i).eq.0) then
         ien2(i) = 0
      else
         ien2(i) = ien1(i) + inc                         
      endif
  100 continue                                                          
c                                                                       
      return                                                            
      end                                                               
c******************************************************************************
      subroutine genside(idside,nside,nen)                            
c                                                                       
c.... program to read and generate element node and material numbers    
c                                                                       
c         idside(nside,nnods) = element sides                        
c                                                                       
      dimension idside(nside,*)                            
c
c.....define idside
c
           if(nen.eq.4) then
              idside(1,1) = 1
              idside(1,2) = 2
c
              idside(2,1) = 2
              idside(2,2) = 3
c
              idside(3,1) = 3
              idside(3,2) = 4
c
              idside(4,1) = 4
              idside(4,2) = 1
           end if
c
            if(nen.eq.9) then
              idside(1,1) = 1
              idside(1,2) = 2
              idside(1,3) = 5
c
              idside(2,1) = 2
              idside(2,2) = 3
              idside(2,3) = 6
c
              idside(3,1) = 3
              idside(3,2) = 4
              idside(3,3) = 7
c
              idside(4,1) = 4
              idside(4,2) = 1
              idside(4,3) = 8
             end if
c
           if(nen.eq.16) then
c
c  lado 1
c
              idside(1,1) = 1
              idside(1,2) = 2
              idside(1,3) = 5
              idside(1,4) = 6
c
c   lado 2
c
              idside(2,1) = 2
              idside(2,2) = 3
              idside(2,3) = 7
              idside(2,4) = 8
c
c   lado 3
c             
              idside(3,1) = 3
              idside(3,2) = 4
              idside(3,3) = 9
              idside(3,4) = 10
c
c   lado 4
c              
              idside(4,1) = 4
              idside(4,2) = 1
              idside(4,3) = 11
              idside(4,4) = 12
c
c
           end if
c
           if(nen.eq.25) then
c
c     lado 1
c
              idside(1,1) = 1
              idside(1,2) = 2
              idside(1,3) = 5
              idside(1,4) = 6
              idside(1,5) = 7
c
c     lado 2
c
              idside(2,1) = 2
              idside(2,2) = 3
              idside(2,3) = 8
              idside(2,4) = 9
              idside(2,5) = 10
c
c     lado 3
c
              idside(3,1) = 3
              idside(3,2) = 4
              idside(3,3) = 11
              idside(3,4) = 12
              idside(3,5) = 13

c
c     lado 4
c
              idside(4,1) = 4
              idside(4,2) = 1
              idside(4,3) = 14
              idside(4,4) = 15
              idside(4,5) = 16

c
           end if
c
c
           if(nen.eq.36) then
c
c     lado 1
c
              idside(1,1) = 1
              idside(1,2) = 2
              idside(1,3) = 5
              idside(1,4) = 6
              idside(1,5) = 7
              idside(1,6) = 8
c
c     lado 2
c
              idside(2,1) = 2
              idside(2,2) = 3
              idside(2,3) = 9
              idside(2,4) = 10
              idside(2,5) = 11
              idside(2,6) = 12
c
c     lado 3
c
              idside(3,1) = 3
              idside(3,2) = 4
              idside(3,3) = 13
              idside(3,4) = 14
              idside(3,5) = 15
              idside(3,6) = 16

c
c     lado 4
c
              idside(4,1) = 4
              idside(4,2) = 1
              idside(4,3) = 17
              idside(4,4) = 18
              idside(4,5) = 19
              idside(4,6) = 20
c
           end if
c
c
           if(nen.eq.49) then
c
c     lado 1
c
              idside(1,1) = 1
              idside(1,2) = 2
              idside(1,3) = 5
              idside(1,4) = 6
              idside(1,5) = 7
              idside(1,6) = 8
              idside(1,7) = 9
c
c     lado 2
c
              idside(2,1) = 2
              idside(2,2) = 3
              idside(2,3) = 10
              idside(2,4) = 11
              idside(2,5) = 12
              idside(2,6) = 13
              idside(2,7) = 14
c
c     lado 3
c
              idside(3,1) = 3
              idside(3,2) = 4
              idside(3,3) = 15
              idside(3,4) = 16
              idside(3,5) = 17
              idside(3,6) = 18
              idside(3,7) = 19

c
c     lado 4
c
              idside(4,1) = 4
              idside(4,2) = 1
              idside(4,3) = 20
              idside(4,4) = 21
              idside(4,5) = 22
              idside(4,6) = 23
              idside(4,7) = 24
c
           end if
c
c
           if(nen.eq.64) then
c
c     lado 1
c
              idside(1,1) = 1
              idside(1,2) = 2
              idside(1,3) = 5
              idside(1,4) = 6
              idside(1,5) = 7
              idside(1,6) = 8
              idside(1,7) = 9
              idside(1,8) = 10
c
c     lado 2
c
              idside(2,1) = 2
              idside(2,2) = 3
              idside(2,3) = 11
              idside(2,4) = 12
              idside(2,5) = 13
              idside(2,6) = 14
              idside(2,7) = 15
              idside(2,8) = 16
c
c     lado 3
c
              idside(3,1) = 3
              idside(3,2) = 4
              idside(3,3) = 17
              idside(3,4) = 18
              idside(3,5) = 19
              idside(3,6) = 20
              idside(3,7) = 21
              idside(3,8) = 22

c
c     lado 4
c
              idside(4,1) = 4
              idside(4,2) = 1
              idside(4,3) = 23
              idside(4,4) = 24
              idside(4,5) = 25
              idside(4,6) = 26
              idside(4,7) = 27
              idside(4,8) = 28
c
           end if
c
c       TRIANGULOS
c
          if(nen.eq.3) then
             idside(1,1) = 1
             idside(1,2) = 2
             idside(2,1) = 2
             idside(2,2) = 3
             idside(3,1) = 3
             idside(3,2) = 1
          end if
c
          if(nen.eq.6) then
             idside(1,1) = 1
             idside(1,2) = 2
             idside(1,3) = 4
            
             idside(2,1) = 2
             idside(2,2) = 3
             idside(2,3) = 5
            
             idside(3,1) = 3
             idside(3,2) = 1
             idside(3,3) = 6
          end if
c
          if(nen.eq.10) then
c         
             idside(1,1) = 1
             idside(1,2) = 2
             idside(1,3) = 4
             idside(1,4) = 5
c          
             idside(2,1) = 2
             idside(2,2) = 3
             idside(2,3) = 6
             idside(2,4) = 7
c             
             idside(3,1) = 3
             idside(3,2) = 1
             idside(3,3) = 8
             idside(3,4) = 9
c
c
          end if
c
          if(nen.eq.15) then
c         
             idside(1,1) = 1
             idside(1,2) = 2
             idside(1,3) = 4
             idside(1,4) = 5
             idside(1,5) = 6
c          
             idside(2,1) = 2
             idside(2,2) = 3
             idside(2,3) = 7
             idside(2,4) = 8
             idside(2,5) = 9
c             
             idside(3,1) = 3
             idside(3,2) = 1
             idside(3,3) = 10
             idside(3,4) = 11
             idside(3,5) = 12
c
          end if
c          
          if(nen.eq.21) then
c         
             idside(1,1) = 1
             idside(1,2) = 2
             idside(1,3) = 4
             idside(1,4) = 5
             idside(1,5) = 6
             idside(1,6) = 7
c          
             idside(2,1) = 2
             idside(2,2) = 3
             idside(2,3) = 8
             idside(2,4) = 9
             idside(2,5) = 10
             idside(2,6) = 11
c             
             idside(3,1) = 3
             idside(3,2) = 1
             idside(3,3) = 12
             idside(3,4) = 13
             idside(3,5) = 14
             idside(3,6) = 15
c
          end if
          
       return 
c                                                                       
      end                
c******************************************************************************
      subroutine genelad(lado,nside)                                     
c                                                                       
c.... program to read and generate element node and material numbers    
c                                                                       
c         lado(nside,numel) = element node numbers                         
c         mat(numel) e,    = element material numbers                     
c         nen            = number of element nodes (le.27)              
c         n              = element number                               
c         ng             = generation parameter                         
c         nel(i)         = number of elements in direction i            
c         incel(i)       = element number increment for direction i     
c         inc(i)         = node number increment for direction i        
c                                                                       
      dimension lado(nside,*),itemp(27)                             
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
      common /genelc/ n,nel(3),incel(3),inc(3)                          
c
  100 continue
      read(iin,1000) n,(itemp(i),i=1,nside),ng  
      if (n.eq.0) return  
cc      write(*,*) n, itemp(1),itemp(2),itemp(3),itemp(4),ng
c
	do no =1,nside
	  lado(no,n) = itemp(no) 
	end do
c
      if (ng.ne.0) then
c                                                                       
c....... generate data                                                     
c                                                                       
         read(iin,1000) (nel(i),incel(i),inc(i),i=1,3)  
c         write(*,*) nel(1),incel(1),inc(1),nel(2),incel(2),inc(2)
         call genelad1(lado,nside)                                          
      endif
c
      go to 100                                                         
c
c                                                                       
 1000 format(20i10)                                             
c                                                                       
      end                                                               
c******************************************************************************
      subroutine genelad1(ipar,nodsp)                                    
c                                                                       
c.... program to generate element node and material numbers             
c                                                                       
      dimension ipar(nodsp,*)                                    
      common /genelc/ n,nel(3),incel(3),inc(3)                          
c
c                                                                       
c.... set defaults                                                      
c                                                                       
      call geneld                                                       
c                                                                       
c.... generation algorithm                                              
c                                                                       
      ie = n                                                            
      je = n                                                            
      ke = n                                                            
c                                                                      
      ii = nel(1)                                                       
      jj = nel(2)                                                       
      kk = nel(3)                                                       
c   
c
      do 300 k=1,kk                                                     
c
      do 200 j=1,jj                                                     
c
      do 100 i=1,ii                                                     
c                                                                       
      if (i.ne.ii) then
         le = ie                                                           
         ie = le + incel(1)                                                
       call genelip(ipar(1,ie),ipar(1,le),inc(1),nodsp)                       
      endif
  100 continue                                                          
c                                                                       
      if (j.ne.jj) then
         le = je                                                           
         je = le + incel(2)                                                
       call genelip(ipar(1,je),ipar(1,le),inc(2),nodsp)                       
         ie = je                                                           
      endif
  200 continue                                                          
c                                                                       
      if (k.ne.kk) then
         le = ke                                                           
         ke = le + incel(3)                                                
       call genelip(ipar(1,ke),ipar(1,le),inc(3),nodsp)                       
         ie = ke                                                           
      endif
  300 continue                                                          
c                                                                       
      return                                                            
      end                                                               
c******************************************************************************
      subroutine genelip(ien2,ien1,inc,nodsp)                              
c                                                                       
c.... program to increment element node numbers                         
c                                                                       
      dimension ien1(*),ien2(*)                                         
c
c                                                                       
      do 100 i=1,nodsp                                                    
      if (ien1(i).eq.0) then
         ien2(i) = 0
      else
         ien2(i) = ien1(i) + inc                         
      endif
  100 continue                                                          
c                                                                       
      return                                                            
      end                                                               
c******************************************************************************
      subroutine genelpar(ipar,ien,lado,idside,
     &              nen,nside,nodsp,numel,npars)      
cc      call genelpar(ipar,ien,lado,idside,
cc     &              nen,nside,nodsp,numel,npars)
c                                                                       
c.... gera numeração dos parametros dos multiplicadores poe elemento   
c                                                                       
c         ipar(npars*nside,numel) = element parameter numbers 
c         nside          = number of element sides                        
c         nodsp          = number of element parameters (le.27)              
c         n              = element number                               
c         ng             = generation parameter                         
c         nel(i)         = number of elements in direction i            
c         incel(i)       = element number increment for direction i     
c         inc(i)         = node number increment for direction i        
c                                                                       
      dimension lado(nside,*),ipar(nodsp,*) 
      dimension ien(nen,*),idside(nside,*)                           
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
c   	
      do 2002 nl=1,numel     
	do 2001 l=1,nside
c
	  la = lado(l,nl)
	  lada = npars*(la-1)

       do np=1,npars
	     nsp = (l-1)*npars + np
	     ipar(nsp,nl) = lada + np
	 end do
c
 2001  continue
 2002  continue 
       return
c                                                                      
c                                                                       
      end                                                               
      
c***********************************************************************      
c     n o d a l  b o d y  f o r c e s      
c***********************************************************************
      subroutine input(f,ndof,numnp,j,nlvect,iprtin)
c
c.... program to read, generate and write nodal input data
c
c        f(ndof,numnp,nlvect) = prescribed forces/kinematic data (j=0)
c                             = nodal body forces(j=1)
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      logical lzero
      dimension f(ndof,numnp,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
      call clear(f,nlvect*numnp*ndof)
c
      do 100 nlv=1,nlvect
      call genfl(f(1,1,nlv),ndof)
      call ztest(f(1,1,nlv),ndof*numnp,lzero)
c
      if (iprtin.eq.0) then
c
         if (lzero) then
            if (j.eq.0) write(iecho,1000) nlv
            if (j.eq.1) write(iecho,2000)
         else
            if (j.eq.0) call printf(f,ndof,numnp,nlv)
c
            if (j.eq.1) 
     &      call printd(' n o d a l  b o d y  f o r c e s  ',
     &                  f,ndof,numnp,iecho)
c
         endif
      endif
c
  100 continue
c
      return
 1000 format(/////,' there are no nonzero prescribed forces and ',
     &    'kinematic boundary conditions for load vector number ',i10)
 2000 format(/////,' there are no nonzero nodal body forces')
      end
c**** new **********************************************************************
      subroutine load(id,f,brhs,ndof,numnp,nlvect)
c
c.... program to accumulate nodal forces and transfer into
c        right-hand-side vector
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension id(ndof,*),f(ndof,numnp,*),brhs(*)
c
      do 300 i=1,ndof
c
      do 200 j=1,numnp
      k = id(i,j)
      if (k.gt.0) then
c
         do 100 nlv=1,nlvect
         brhs(k) = brhs(k) + f(i,j,nlv)
  100    continue
c
      endif
c
  200 continue
c
  300 continue
c
      return
      end
c**** new **********************************************************************
      subroutine ftod(id,d,f,ndof,numnp,nlvect)
c
c.... program to compute displacement boundary conditions
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension id(ndof,*),d(ndof,*),f(ndof,numnp,*)
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      do 300 i=1,ndof
c
            do 200 j=1,numnp
c
            k = id(i,j)
            if (k.gt.0) go to 200
            val = zero
                  do 100 lv=1,nlvect
                  val = val + f(i,j,lv)
100               continue
c
            d(i,j) = val
c
  200       continue
c
 300  continue
      return
      end
      
c***********************************************************************
c          shape functions
c***********************************************************************
      subroutine oneshl(shl,w,nint,nen)
c
c.... program to calculate integration-rule weights, shape functions 
c        and local derivatives for a one-dimensional element
c
c                 r = local element coordinate ("xi")
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = shape function
c              w(l) = integration-rule weight
c                 i = local node number 
c                 l = integration-point number
c              nint = number of integration points, eq. 1, 2, 3, 4, 6, 7 ,8
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension shl(3,nen,*),w(*),ra(10),xa(10)
      data   zero,pt1667,pt25,pt5 
     &      /0.0d0,0.1666666666666667d0,0.25d0,0.5d0/, 
     &       one,two,three,four,five,six 
     &      /1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0/ 

      data  five9/0.5555555555555555d0/,eight9/0.8888888888888888d0/
c
      if (nint.eq.1) then
         w(1)  = two
         ra(1) = zero
      endif
c
      if (nint.eq.2) then
         w(1) = one
         w(2) = one
         ra(1)=-.577350269189626
         ra(2)= .577350269189625
      endif
c
      if (nint.eq.3) then
         w(1) = five9
         w(2) = eight9
         w(3) = five9
         ra(1)=-.774596669241483
         ra(2)= zero
         ra(3)= .774596669241483
      endif
c
      if (nint.eq.4) then
         w(1) = .347854845137454
         w(2) = .652145154862546
         w(3) = .652145154862546
         w(4) = .347854845137454
         ra(1)=-.861136311594053
         ra(2)=-.339981043584856
         ra(3)= .339981043584856
         ra(4)= .861136311594053
      endif
c
       if(nint.eq.5) then
        w(1) = .236926885056189
        w(2) = .478628670499366
        w(3) = .568888888888888
        w(4) = .478628670499366
        w(5) = .236926885056189
        ra(1)=-.906179845938664
        ra(2)=-.538469310105683
        ra(3)= zero
        ra(4)= .538469310105683
        ra(5)= .906179845938664
       endif
c
       if(nint.eq.6) then
         w(1) = .171324492397170
         w(2) = .360761573048139
         w(3) = .467913934572691
         w(4) = .467913934572691
         w(5) = .360761573048139
         w(6) = .171324492397170
c         
         ra(1)=-.932469514203152
         ra(2)=-.661209386466265
         ra(3)=-.238619186083197
         ra(4)= .238619186083197
         ra(5)= .661209386466365
         ra(6)= .932469514203152
        endif
c
       if(nint.eq.7) then
         w(1) = .129484966168870
         w(2) = .279705391489277 
         w(3) = .381830050505119
         w(4) = .417959183673469
         w(5) = .381830050505119
         w(6) = .279705391489277
         w(7) = .129484966168870
c         
         ra(1)=-.949107912342759
         ra(2)=-.741531185599394
         ra(3)=-.405845151377397
         ra(4)= zero
         ra(5)= .405845151377397
         ra(6)= .741531185599394
         ra(7)= .949107912342759
        endif
c
       if(nint.eq.8) then
         w(1) = .101228536290376
         w(2) = .222381034453374
         w(3) = .313706645877887
         w(4) = .362683783378362
         w(5) = .362683783378362
         w(6) = .313706645877887
         w(7) = .222381034453374
         w(8) = .101228536290376
c
         ra(1)=-.960289856497536
         ra(2)=-.796666477413627
         ra(3)=-.525532409916329
         ra(4)=-.183434642495650
         ra(5)= .183434642495650
         ra(6)= .525532409916329
         ra(7)= .796666477413627
         ra(8)= .960289856497536
        endif

      if (nen.eq.1) xa(1) = zero
c
      if (nen.eq.2) then
         xa(1) = -one
         xa(2) =  one
      endif
c 
      if(nen.eq.3) then
         xa(1)= -one
         xa(2)= one
         xa(3)= zero
      endif
c
      if (nen.eq.4) then
         xa(1) = -one
         xa(2) = one
         xa(3) = -.333333333333333
         xa(4) =  .333333333333333
         endif
c
       if(nen.eq.5) then
         xa(1)= -one 
         xa(2)=  one 
         xa(3)= -pt5
         xa(4)= zero
         xa(5)= pt5
       endif
c
        if(nen.eq.6) then
         xa(1) = -one
         xa(2) =  one
         xa(3) = -.600000000000000
         xa(4) = -.200000000000000
         xa(5) =  .200000000000000
         xa(6) =  .600000000000000
        endif
c
        if(nen.eq.7) then
         xa(1) = -one
         xa(2) =  one
         xa(3) = -.666666666666666
         xa(4) = -.333333333333333
         xa(5) = zero
         xa(6) =  .333333333333333
         xa(7) =  .666666666666666
        endif

        if(nen.eq.8) then
         xa(1) = -one
         xa(2) =  one
         xa(3) = -0.71428571428571
         xa(4) = -0.42857142857143
         xa(5) = -0.14285714285714
         xa(6) = 0.14285714285714
         xa(7) = 0.42857142857143
         xa(8) = 0.71428571428571
        endif
c
      do 100 l = 1, nint
         r = ra(l)
c
        if(nen.eq.1) then
        shl(1,1,l) = zero
        shl(2,1,l) = one
        go to 100
        endif
c
        do 50 i = 1, nen
         aa = one
         bb = one
         aax = zero
         do 40 j =1, nen
          daj = one
          if (i .ne. j)then
           aa = aa * ( r - xa(j))
           bb = bb * ( xa(i) - xa(j))
           do 30 k = 1, nen
            if(k .ne. i .and. k .ne. j) then
              daj = daj * ( r - xa(k))
            endif  
   30      continue
           aax =aax + daj
          endif
   40    continue
        shl(2,i,l) = aa/bb
        shl(1,i,l) = aax/bb
   50  continue
c
  100  continue      
       return
       end
c******* nova subrotina ***************************************************************
c                                              
      subroutine shlt(shl,w,nint,nen) 
c
c     Objetivo: calcular pesos, funcoes de interpolacao e derivadas locais
c               para elementos triangulares
c
c----------------------------------------------------------------------
c
c 
c.... program to calculate integration-rule weights, shape functions 
c        and local derivatives for a triangular element 
c 
c        c1, c2, c3 = local element coordinates ("l1", "l2", "l3".) 
c        shl(j,i,l) = local ("j") derivative of shape function 
c        shl(3,i,l) = local  shape function 
c              w(l) = integration-rule weight 
c                 i = local node number 
c                 l = integration point number 
c              nint = number of integration points,
c 
      implicit real*8 (a-h,o-z) 
c 
c.... remove above card for single precision operation 
c 
      dimension shl(3,nen,*),w(*),cl1(25),cl2(25),cl3(25)
      data   zero,pt1667,pt25,pt5 
     &      /0.0d0,0.1666666666666667d0,0.25d0,0.5d0/, 
     &       one,two,three,four,five,six 
     &      /1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0/ 
      data r1/0.33333333333333333333d00/,w1/1.d00/, 
     &     r2/0.5d00                   /,w2/0.3333333333333333333d00/, 
     &     r3a/0.3333333333333333333d00/,w3a/-0.5625d00/, 
     &     r3b1/0.6d00                 /,w3b/0.520833333333333333d00/, 
     &     r3b2/0.2d00                 /, 
     &     r7a/0.3333333333333333333d00/,w7a/0.225d00/, 
     &     r7b/0.0597158717d00         /,w7b/0.1323941527d00/, 
     &     r7c/0.4701420641d00         /, 
     &     r7d/0.7974269853d00         /,w7d/0.1259391805d00/, 
     &     r7e/0.1012865073d00         /
      data ri1/0.1666666666666666666d00/,ri2/0.666666666666666666d00/,
     &     ri3/0.1666666666666666666d00/
      data r2a/0.666666666666667d00    /,r2b/ 0.166666666666667d00/,
     &     w6a/0.109951743655322d00    /,w6b/0.223381589678011d00/,
     &     r6a/0.816847572980459d00    /,r6c/0.108103018168070d00/,
     &     r6b/0.091576213509771d00    /,r6d/0.445948490915965d00/,
     &     w12a/0.116786275726379d00   /,
     &     r12a1/0.501426509658179d00  /,r12b1/0.249286745170910d00/,
     &     w12b/0.050844906370207d00   /,
     &     r12a2/0.873821971016996d00  /,r12b2/0.063089014491502d00/,
     &     w12c/0.082851075618374d00   /,
     &     r12a3/0.636502499121399d00  /,r12b3/0.310352451033785d00/,
     &     r12c3/0.053145049844816d00  /,
c     
     &     w13a/-0.149570044467670d00  /,w13b/0.175615257433204d00/,
     &     w13c/0.053347235608839d00   /,w13d/0.077113760890257d00 /,
     &     r13a/0.333333333333333d00/,
     &     r13a1/0.479308067841923d00  /,r13b1/0.260345966079038d00/,
     &     r13a2/0.869739794195568d00  /,r13b2/0.065130102902216d00/,
     &     r13a3/0.638444188569809d00  /,r13b3/0.312865496004875d00/,
     &     r13c3/0.048690315425316d00  /
      data w16a/0.144315607677787d00   /,r16a/0.333333333333333d00/,
     &     w16b/0.095091634267285d00   /,
     &     r16b1/0.081414823414554d00  /,r16b2/0.459292588292723d00/,
     &     w16c/0.103217370534718d00   /,
     &     r16c1/0.658861384496480d00  /,r16c2/0.170569307751760d00/,
     &     w16d/0.032458497623198d00   /,
     &     r16d1/0.898905543365938d00  /,r16d2/0.050547228317031d00/,
     &     w16e/0.027230314174435d00   /,r16e1/0.728492392955404d00/,
     &     r16e2/0.263112829634638d00  /,r16e3/0.008394777409958d00/,
c     
     &     w19a/0.097135796282799d00   /,r19a/0.333333333333333d00/, 
     &     w19b/0.031334700227139d00   /,
     &     r19b1/0.020634961602525d00  /,r19b2/0.489682519198738d00/, 
     &     w19c/0.077827541004774d00   /,
     &     r19c1/0.125820817014127d00  /,r19c2/0.437089591492937d00/, 
     &     w19d/0.079647738927210d00   /,
     &     r19d1/0.623592928761935d00  /,r19d2/0.188203535619033d00/, 
     &     w19e/0.025577675658698d00   /,
     &     r19e1/0.910540973211095d00  /,r19e2/0.044729513394453d00/, 
     &     w19f/0.043283539377289d00   /,r19f1/0.741198598784498d00/,
     &     r19f2/0.221962989160766d00  /,r19f3/0.036838412054736d00/,
c
     &     w25a/0.090817990382754d00   /,r25a/0.333333333333333d00/, 
     &     w25b/0.036725957756467d00   /,
     &     r25b1/0.028844733232685d00  /,r25b2/0.485577633383657d00/, 
     &     w25c/0.045321059435528d00   /,
     &     r25c1/0.781036849029926d00  /,r25c2/0.109481575485037d00/, 
     &     w25d/0.072757916845420d00   /,r25d1/0.550352941820999d00/,
     &     r25d2/0.307939838764121d00  /,r25d3/0.141707219414880d00/,
     &     w25e/0.028327242531057d00   /,r25e1/0.728323904597411d00/,
     &     r25e2/0.246672560639903d00  /,r25e3/0.025003534762686d00/,
     &     w25f/0.009421666963733d00   /,r25f1/0.923655933587500d00/,
     &     r25f2/0.066803251012200d00  /,r25f3/0.009540815400299d00/

      pt3s2 = three/two 
      pt9s2 = 9.d0/two
      pt27s2= 27.d0/two
      pt2s3 = two/three
      pt1s6 = one/six
      pt32s3= 32.d0/three
      pt8s3 = 8.d0/three
c 
      if (nint.eq.1) then 
            w(1)=w1/two 
            cl1(1)=r1 
            cl2(1)=r1 
            cl3(1)=one-r1-r1 
      end if 
c 
      if(nint.eq.3) then 
        do 31 i=1,3 
          w(i)=w2/two 
  31    continue 
        cl1(1)=r2 
        cl2(1)=r2 
        cl3(1)=zero 
        cl1(2)=zero 
        cl2(2)=r2 
        cl3(2)=r2 
        cl1(3)=r2 
        cl2(3)=zero 
        cl3(3)=r2 
! c  pontos montados nas linhas do meio (artigo Dunavant)
! c           ponto a
!             cl1(1)=r2a
!             cl2(1)=r2b
!             cl3(1)=r2b
! c           ponto b
!             cl1(2)=r2b
!             cl2(2)=r2a
!             cl3(2)=r2b
! c           ponto c
!             cl1(3)=r2b
!             cl2(3)=r2b
!             cl3(3)=r2a
! 
      end if 
c 
      if(nint.eq.4) then 
        w(1)= w3a/two 
        do 41 i=2,4 
          w(i)=w3b/two 
  41    continue 
        cl1(1)=r3a 
        cl2(1)=r3a 
        cl3(1)=r3a
        
        cl1(2)=r3b1 
        cl2(2)=r3b2 
        cl3(2)=r3b2 
        
        cl1(3)=r3b2 
        cl2(3)=r3b1 
        cl3(3)=r3b2 
        
        cl1(4)=r3b2 
        cl2(4)=r3b2 
        cl3(4)=r3b1 
      end if 
c
      if(nint.eq.6) then
        do 61 i=1,3
          w(i)=w6a/two
  61    continue
        do 62 i=4,6
          w(i)=w6b/two
  62    continue
        do 63 i=1,3
          cl1(i)=r6b
          cl2(i)=r6b
          cl3(i)=r6b
  63    continue
          cl1(1)=r6a
          cl2(2)=r6a
          cl3(3)=r6a
        do 64 i=4,6
          cl1(i)=r6d
          cl2(i)=r6d
          cl3(i)=r6d
  64    continue
          cl1(4)=r6c
          cl2(5)=r6c
          cl3(6)=r6c
      end if
c
      if(nint.eq.7) then 
        w(1)= w7a/two 
        do 71 i=2,4 
          w(i)=w7b/two 
  71    continue 
        do 72 i=5,7 
          w(i)=w7d/two 
  72    continue 
          cl1(1)=r7a 
          cl2(1)=r7a 
          cl3(1)=r7a 
        do 73 i=2,4 
          cl1(i)=r7c 
          cl2(i)=r7c 
          cl3(i)=r7c 
  73    continue 
          cl1(2)=r7b 
          cl2(3)=r7b 
          cl3(4)=r7b 
        do 74 i=5,7 
          cl1(i)=r7e 
          cl2(i)=r7e 
          cl3(i)=r7e 
  74    continue 
          cl1(5)=r7d 
          cl2(6)=r7d 
          cl3(7)=r7d 
      end if 
c
      if(nint.eq.12) then
          do 110 i=1,3
            w(i)=w12a/two
  110     continue
          do 120 i=4,6
            w(i)=w12b/two
  120     continue
          do 130 i=7,12
            w(i)=w12c/two
  130     continue
          do 140 i=1,3
            cl1(i)=r12b1
            cl2(i)=r12b1
            cl3(i)=r12b1
  140     continue
          cl1(1)=r12a1
          cl2(2)=r12a1
          cl3(3)=r12a1
          do 150 i=4,6
            cl1(i)=r12b2
            cl2(i)=r12b2
            cl3(i)=r12b2
  150     continue
          cl1(4)=r12a2
          cl2(5)=r12a2
          cl3(6)=r12a2
          
          cl1(7)=r12a3
          cl2(7)=r12b3
          cl3(7)=r12c3
          
          cl1(8)=r12c3
          cl2(8)=r12a3
          cl3(8)=r12b3
          
          cl1(9)=r12b3
          cl2(9)=r12c3
          cl3(9)=r12a3
          
          cl1(10)=r12a3
          cl2(10)=r12c3
          cl3(10)=r12b3
          
          cl1(11)=r12b3
          cl2(11)=r12a3
          cl3(11)=r12c3
          
          cl1(12)=r12c3
          cl2(12)=r12b3
          cl3(12)=r12a3
      end if
cc
      if (nint.eq.13) then         
        w(1)= w13a/two 
        do 131 i=2,4 
          w(i)=w13b/two 
  131   continue 
        do 132 i=5,7 
          w(i)=w13c/two 
  132   continue 
        do 133 i=8,13 
          w(i)=w13d/two 
  133   continue 
          cl1(1)=r13a 
          cl2(1)=r13a 
          cl3(1)=r13a 
        do 134 i=2,4 
          cl1(i)=r13b1 
          cl2(i)=r13b1
          cl3(i)=r13b1 
  134   continue 
          cl1(2)=r13a1 
          cl2(3)=r13a1 
          cl3(4)=r13a1 
        do 135 i=5,7 
          cl1(i)=r13b2 
          cl2(i)=r13b2 
          cl3(i)=r13b2 
  135   continue 
          cl1(5)=r13a2 
          cl2(6)=r13a2 
          cl3(7)=r13a2 
          
          cl1(8)=r13a3 
          cl2(8)=r13b3 
          cl3(8)=r13c3 
          
          cl1(9)=r13c3 
          cl2(9)=r13a3 
          cl3(9)=r13b3 
          
          cl1(10)=r13b3 
          cl2(10)=r13c3 
          cl3(10)=r13a3 
 
          cl1(11)=r13a3 
          cl2(11)=r13c3 
          cl3(11)=r13b3 

          cl1(12)=r13b3 
          cl2(12)=r13a3 
          cl3(12)=r13c3 
          
          cl1(13)=r13c3 
          cl2(13)=r13b3 
          cl3(13)=r13a3 
      endif 
c
      if (nint.eq.16) then         
        w(1)= w16a/two 
        do 161 i=2,4 
          w(i)=w16b/two 
  161   continue 
        do 162 i=5,7 
          w(i)=w16c/two 
  162   continue 
        do 163 i=8,10 
          w(i)=w16d/two 
  163   continue 
        do 164 i=11,16 
          w(i)=w16e/two 
  164   continue 
          cl1(1)=r16a 
          cl2(1)=r16a 
          cl3(1)=r16a 
        do 165 i=2,4 
          cl1(i)=r16b2 
          cl2(i)=r16b2
          cl3(i)=r16b2 
  165   continue 
          cl1(2)=r16b1 
          cl2(3)=r16b1 
          cl3(4)=r16b1 
        do 166 i=5,7 
          cl1(i)=r16c2 
          cl2(i)=r16c2 
          cl3(i)=r16c2 
  166   continue 
          cl1(5)=r16c1 
          cl2(6)=r16c1 
          cl3(7)=r16c1 
        do 167 i=8,10 
          cl1(i)=r16d2 
          cl2(i)=r16d2 
          cl3(i)=r16d2 
  167   continue 
          cl1(8)=r16d1 
          cl2(9)=r16d1 
          cl3(10)=r16d1 
          
          cl1(11)=r16e1 
          cl2(11)=r16e2 
          cl3(11)=r16e3 
          
          cl1(12)=r16e3 
          cl2(12)=r16e1 
          cl3(12)=r16e2 
          
          cl1(13)=r16e2 
          cl2(13)=r16e3 
          cl3(13)=r16e1 
 
          cl1(14)=r16e1 
          cl2(14)=r16e3 
          cl3(14)=r16e2 

          cl1(15)=r16e2 
          cl2(15)=r16e1 
          cl3(15)=r16e3 
          
          cl1(16)=r16e3 
          cl2(16)=r16e2 
          cl3(16)=r16e1 
      endif 
c
      if (nint.eq.19) then         
        w(1)= w19a/two 
        do 191 i=2,4 
          w(i)=w19b/two 
  191   continue 
        do 192 i=5,7 
          w(i)=w19c/two 
  192   continue 
        do 193 i=8,10 
          w(i)=w19d/two 
  193   continue 
        do 194 i=11,13 
          w(i)=w19e/two 
  194   continue 
        do 195 i=14,19 
          w(i)=w19f/two 
  195   continue 
          cl1(1)=r19a 
          cl2(1)=r19a 
          cl3(1)=r19a 
        do 196 i=2,4 
          cl1(i)=r19b2 
          cl2(i)=r19b2
          cl3(i)=r19b2 
  196   continue 
          cl1(2)=r19b1 
          cl2(3)=r19b1 
          cl3(4)=r19b1 
        do 197 i=5,7 
          cl1(i)=r19c2 
          cl2(i)=r19c2 
          cl3(i)=r19c2 
  197   continue 
          cl1(5)=r19c1 
          cl2(6)=r19c1 
          cl3(7)=r19c1 
        do 198 i=8,10 
          cl1(i)=r19d2 
          cl2(i)=r19d2 
          cl3(i)=r19d2 
  198   continue 
          cl1(8)=r19d1 
          cl2(9)=r19d1 
          cl3(10)=r19d1 
        do 199 i=11,13 
          cl1(i)=r19e2 
          cl2(i)=r19e2 
          cl3(i)=r19e2 
  199   continue 
          cl1(11)=r19e1 
          cl2(12)=r19e1 
          cl3(13)=r19e1 
          
          cl1(14)=r19f1 
          cl2(14)=r19f2 
          cl3(14)=r19f3 
          
          cl1(15)=r19f3 
          cl2(15)=r19f1 
          cl3(15)=r19f2 
          
          cl1(16)=r19f2 
          cl2(16)=r19f3 
          cl3(16)=r19f1 
 
          cl1(17)=r19f1 
          cl2(17)=r19f3 
          cl3(17)=r19f2 

          cl1(18)=r19f2 
          cl2(18)=r19f1 
          cl3(18)=r19f3 
          
          cl1(19)=r19f3 
          cl2(19)=r19f2 
          cl3(19)=r19f1 
      endif 
c
      if (nint.eq.25) then         
        w(1)= w25a/two 
        do 251 i=2,4 
          w(i)=w25b/two 
  251   continue 
        do 252 i=5,7 
          w(i)=w25c/two 
  252   continue 
        do 253 i=8,13 
          w(i)=w25d/two 
  253   continue 
        do 254 i=14,19 
          w(i)=w25e/two 
  254   continue 
        do 255 i=20,25 
          w(i)=w25f/two 
  255   continue 
          cl1(1)=r25a 
          cl2(1)=r25a 
          cl3(1)=r25a 
        do 256 i=2,4 
          cl1(i)=r25b2 
          cl2(i)=r25b2
          cl3(i)=r25b2 
  256   continue 
          cl1(2)=r25b1 
          cl2(3)=r25b1 
          cl3(4)=r25b1 
        do 257 i=5,7 
          cl1(i)=r25c2 
          cl2(i)=r25c2 
          cl3(i)=r25c2 
  257   continue 
          cl1(5)=r25c1 
          cl2(6)=r25c1 
          cl3(7)=r25c1 
          
          cl1(8)=r25d1 
          cl2(8)=r25d2 
          cl3(8)=r25d3 
          
          cl1(9)=r25d3 
          cl2(9)=r25d1 
          cl3(9)=r25d2 
          
          cl1(10)=r25d2 
          cl2(10)=r25d3 
          cl3(10)=r25d1 
 
          cl1(11)=r25d1 
          cl2(11)=r25d3 
          cl3(11)=r25d2 

          cl1(12)=r25d2 
          cl2(12)=r25d1 
          cl3(12)=r25d3 
          
          cl1(13)=r25d3 
          cl2(13)=r25d2 
          cl3(13)=r25d1 

          cl1(14)=r25e1 
          cl2(14)=r25e2 
          cl3(14)=r25e3 
          
          cl1(15)=r25e3 
          cl2(15)=r25e1 
          cl3(15)=r25e2 
          
          cl1(16)=r25e2 
          cl2(16)=r25e3 
          cl3(16)=r25e1 
 
          cl1(17)=r25e1 
          cl2(17)=r25e3 
          cl3(17)=r25e2 

          cl1(18)=r25e2 
          cl2(18)=r25e1 
          cl3(18)=r25e3 
          
          cl1(19)=r25e3 
          cl2(19)=r25e2 
          cl3(19)=r25e1 
          
          cl1(20)=r25f1 
          cl2(20)=r25f2 
          cl3(20)=r25f3 
          
          cl1(21)=r25f3 
          cl2(21)=r25f1 
          cl3(21)=r25f2 
          
          cl1(22)=r25f2 
          cl2(22)=r25f3 
          cl3(22)=r25f1 
 
          cl1(23)=r25f1 
          cl2(23)=r25f3 
          cl3(23)=r25f2 

          cl1(24)=r25f2 
          cl2(24)=r25f1 
          cl3(24)=r25f3 
          
          cl1(25)=r25f3 
          cl2(25)=r25f2 
          cl3(25)=r25f1 

      endif 

      do 200 l=1,nint 
c 
            c1 = cl1(l) 
            c2 = cl2(l) 
            c3 = cl3(l) 
!          acrescentei o if (2012-10-08)
!           interpolação linear (p=1 e nen=3)(2012-10-08)
            if(nen.eq.1) then
              shl(1,1,l)= zero
              shl(2,1,l)= zero
              shl(3,1,l)= one
              shl(1,2,l)= zero
              shl(2,2,l)= zero
              shl(3,2,l)= one
              shl(1,3,l)= zero
              shl(2,3,l)= zero
              shl(3,3,l)= one
            end if 
            if(nen.eq.3) then
              shl(1,1,l)= one
              shl(2,1,l)= zero
              shl(3,1,l)= c1
              shl(1,2,l)= zero
              shl(2,2,l)= one
              shl(3,2,l)= c2
              shl(1,3,l)=-one
              shl(2,3,l)=-one
              shl(3,3,l)= c3
            end if 
!           interpolação quadrática (p=2 e nen=6)(2012-10-08)
            if(nen.eq.6) then
              shl(1,1,l)= four*c1-one
              shl(2,1,l)= zero
              shl(3,1,l)= (two*c1 - one)*c1
              shl(1,2,l)= zero
              shl(2,2,l)= four*c2-one
              shl(3,2,l)= (two*c2 - one)*c2
              shl(1,3,l)= one - four*c3
              shl(2,3,l)= one - four*c3
              shl(3,3,l)= (two*c3 - one)*c3
              shl(1,4,l)= four * c2
              shl(2,4,l)= four * c1
              shl(3,4,l)= four * c1 * c2
              shl(1,5,l)=-four * c2
              shl(2,5,l)= four * (c3 - c2)
              shl(3,5,l)= four * c2 * c3
              shl(1,6,l)= four * (c3 - c1)
              shl(2,6,l)=-four * c1
              shl(3,6,l)= four * c3 * c1
            end if
!           interpolação cúbica (p=3 e nen=10)(2012-10-22)
            if(nen.eq.10) then
           
              shl(1,1,l)= pt5*(three*c1-two)*(three*c1-one) +
     &                pt3s2*c1*(three*c1-one)+pt3s2*c1*(three*c1-two)
              shl(2,1,l)= zero
              shl(3,1,l)= pt5*c1*(three*c1 - two)*(three*c1-one)
              
              shl(1,2,l)= zero
              shl(2,2,l)= pt5*(three*c2-two)*(three*c2-one) +
     &                pt3s2*c2*(three*c2-one)+pt3s2*c2*(three*c2-two)
              shl(3,2,l)= pt5*c2*(three*c2 - two)*(three*c2-one)
              
              shl(1,3,l)= -pt5*(three*c3 - two)*(three*c3-one)  
     &          - pt3s2*c3*(three*c3-one) - pt3s2*c3*(three*c3 - two)
              shl(2,3,l)= -pt5*(three*c3 - two)*(three*c3-one)  
     &          - pt3s2*c3*(three*c3-one) - pt3s2*c3*(three*c3 - two)
              shl(3,3,l)= pt5*c3*(three*c3 - two)*(three*c3-one)  
              
              shl(1,4,l)= pt9s2*c2*(three*c1-one) + pt27s2*c1*c2
              shl(2,4,l)= pt9s2*c1*(three*c1-one)
              shl(3,4,l)= pt9s2*c1*c2*(three*c1-one)
              
              shl(1,5,l)= pt9s2*c2*(three*c2-one)
              shl(2,5,l)= pt9s2*c1*(three*c2-one) + pt27s2*c1*c2
              shl(3,5,l)= pt9s2*c1*c2*(three*c2-one)
              
              shl(1,6,l)= -pt9s2*c2*(three*c2-one)
              shl(2,6,l)= pt9s2*c3*(three*c2-one) + pt27s2*c2*c3
     &                  - pt9s2*c2*(three*c2-one)         
              shl(3,6,l)= pt9s2*c2*c3*(three*c2-one)

              shl(1,7,l)= -pt9s2*c2*(three*c3-one) -pt27s2*c2*c3
              shl(2,7,l)= pt9s2*c3*(three*c3-one) - pt27s2*c2*c3
     &                  - pt9s2*c2*(three*c3-one)         
              shl(3,7,l)= pt9s2*c2*c3*(three*c3-one)

              shl(1,8,l)= -pt9s2*c1*(three*c3-one) -pt27s2*c3*c1
     &                  + pt9s2*c3*(three*c3-one)         
              shl(2,8,l)= -pt9s2*c1*(three*c3-one) - pt27s2*c3*c1
              shl(3,8,l)= pt9s2*c3*c1*(three*c3-one)

              shl(1,9,l)= -pt9s2*c1*(three*c1-one) +pt27s2*c3*c1
     &                  + pt9s2*c3*(three*c1-one)         
              shl(2,9,l)= -pt9s2*c1*(three*c1-one)
              shl(3,9,l)= pt9s2*c3*c1*(three*c1-one)

              shl(1,10,l)= 27.d0*c2*c3 - 27.d0*c1*c2
              shl(2,10,l)= 27.d0*c3*c1 - 27.d0*c1*c2
              shl(3,10,l)= 27.d0*c1*c2*c3

            end if
!           interpolação quártica (p=4 e nen=15)(2012-10-23)
            if(nen.eq.15) then
           
              shl(1,1,l)= pt2s3*(four*c1-two)*(four*c1-one)*c1
     &                  + pt2s3*(four*c1-three)*(four*c1-one)*c1
     &                  + pt2s3*(four*c1-three)*(four*c1-two)*c1
     &            + pt1s6*(four*c1-three)*(four*c1-two)*(four*c1-one)
              shl(2,1,l)= zero
              shl(3,1,l)= pt1s6*(four*c1-three)*(four*c1-two)*
     &                  (four*c1-one)*c1
              
              shl(1,2,l)= zero
              shl(2,2,l)= pt2s3*(four*c2-two)*(four*c2-one)*c2
     &                  + pt2s3*(four*c2-three)*(four*c2-one)*c2
     &                  + pt2s3*(four*c2-three)*(four*c2-two)*c2
     &            + pt1s6*(four*c2-three)*(four*c2-two)*(four*c2-one)
              shl(3,2,l)= pt1s6*(four*c2-three)*(four*c2-two)*
     &                  (four*c2-one)*c2
              
              shl(1,3,l)= -pt2s3*(four*c3-two)*(four*c3-one)*c3
     &                  -pt2s3*(four*c3-three)*(four*c3-one)*c3
     &                  -pt2s3*(four*c3-three)*(four*c3-two)*c3
     &             -pt1s6*(four*c3-three)*(four*c3-two)*(four*c3-one)
              shl(2,3,l)= -pt2s3*(four*c3-two)*(four*c3-one)*c3
     &                  -pt2s3*(four*c3-three)*(four*c3-one)*c3
     &                  -pt2s3*(four*c3-three)*(four*c3-two)*c3
     &             -pt1s6*(four*c3-three)*(four*c3-two)*(four*c3-one)
              shl(3,3,l)= pt1s6*(four*c3-three)*(four*c3-two)*
     &                   (four*c3-one)*c3
              
              shl(1,4,l)= pt32s3*c2*(four*c1-one)*c1
     &                  +pt32s3*c2*(four*c1-two)*c1
     &                  +pt8s3*c2*(four*c1-two)*(four*c1-one)
              shl(2,4,l)= pt8s3*(four*c1-two)*(four*c1-one)*c1
              shl(3,4,l)= pt8s3*c2*(four*c1-two)*(four*c1-one)*c1
              
              shl(1,5,l)= 16.d0*c1*(four*c2-one)*c2 
     &                  + four*(four*c1-one)*(four*c2-one)*c2
              shl(2,5,l)= 16.d0*c2*(four*c1-one)*c1
     &                  + four*(four*c1-one)*c1*(four*c2-one)
              shl(3,5,l)= four*c2*(four*c2-one)*(four*c1-one)*c1
              
              shl(1,6,l)= pt8s3*(four*c2-two)*(four*c2-one)*c2
              shl(2,6,l)= pt32s3*c1*(four*c2-one)*c2 
     &                  + pt32s3*c1*(four*c2-two)*c2 
     &                  + pt8s3*c1*(four*c2-two)*(four*c2-one)
              shl(3,6,l)= pt8s3*c1*(four*c2-two)*(four*c2-one)*c2

              shl(1,7,l)= -pt8s3*(four*c2-two)*(four*c2-one)*c2
              shl(2,7,l)= -pt8s3*(four*c2-two)*(four*c2-one)*c2 
     &                  + pt32s3*c3*(four*c2-one)*c2
     &                  + pt32s3*c3*(four*c2-two)*c2
     &                  + pt8s3*c3*(four*c2-two)*(four*c2-one)
              shl(3,7,l)= pt8s3*c3*(four*c2-two)*(four*c2-one)*c2

              shl(1,8,l)= -16.d0*(four*c2-one)*c2*c3
     &                  -four*(four*c2-one)*c2*(four*c3-one)
              shl(2,8,l)= 16.d0*c2*(four*c3-one)*c3
     &                  +four*(four*c2-one)*(four*c3-one)*c3
     &                  -16.d0*(four*c2-one)*c2*c3
     &                  -four*(four*c2-one)*c2*(four*c3-one)
              shl(3,8,l)= four*(four*c2-one)*c2*(four*c3-one)*c3

              shl(1,9,l)= -pt32s3*c2*(four*c3-one)*c3
     &                  -pt32s3*c2*(four*c3-two)*c3
     &                  -pt8s3*c2*(four*c3-two)*(four*c3-one)
              shl(2,9,l)= pt8s3*(four*c3-two)*(four*c3-one)*c3
     &                  -pt32s3*c2*(four*c3-one)*c3
     &                  -pt32s3*c2*(four*c3-two)*c3
     &                  -pt8s3*c2*(four*c3-two)*(four*c3-one)
              shl(3,9,l)= pt8s3*c2*(four*c3-two)*(four*c3-one)*c3

              shl(1,10,l)= pt8s3*(four*c3-two)*(four*c3-one)*c3
     &                   -pt32s3*c1*(four*c3-one)*c3
     &                   -pt32s3*c1*(four*c3-two)*c3
     &                   -pt8s3*c1*(four*c3-two)*(four*c3-one)
              shl(2,10,l)= -pt32s3*c1*(four*c3-one)*c3
     &                   -pt32s3*c1*(four*c3-two)*c3
     &                   -pt8s3*c1*(four*c3-two)*(four*c3-one)
              shl(3,10,l)= pt8s3*c1*(four*c3-two)*(four*c3-one)*c3

              shl(1,11,l)= -16.d0*c3*(four*c1-one)*c1
     &                   -4.d0*(four*c3-one)*(four*c1-one)*c1
     &                   +16.d0*c1*(four*c3-one)*c3
     &                   +4.d0*(four*c3-one)*c3*(four*c1-one)
              shl(2,11,l)= -16.d0*c3*(four*c1-one)*c1
     &                   -4.d0*(four*c3-one)*(four*c1-one)*c1
              shl(3,11,l)= 4.d0*(four*c3-one)*c3*(four*c1-one)*c1
              
              shl(1,12,l)= -pt8s3*(four*c1-two)*(four*c1-one)*c1
     &                   +pt32s3*c3*(four*c1-one)*c1
     &                   +pt32s3*c3*(four*c1-two)*c1
     &                   +pt8s3*c3*(four*c1-two)*(four*c1-one)
              shl(2,12,l)= -pt8s3*(four*c1-two)*(four*c1-one)*c1
              shl(3,12,l)= pt8s3*c3*(four*c1-two)*(four*c1-one)*c1
              
              shl(1,13,l)= 32.d0*c2*c3*(four*c1-one)
     &                   -32.d0*c1*c2*(four*c1-one)+128.d0*c1*c2*c3
              shl(2,13,l)= 32.d0*c1*c3*(four*c1-one)
     &                   -32.d0*c1*c2*(four*c1-one)
              shl(3,13,l)= 32.d0*c1*c2*c3*(four*c1-one)
              
              shl(1,14,l)= 32.d0*(four*c2-one)*c2*c3
     &                   -32.d0*c1*(four*c2-one)*c2
              shl(2,14,l)= 32.d0*c1*c3*(four*c2-one)
     &                   -32.d0*c1*c2*(four*c2-one)+128.d0*c1*c2*c3
              shl(3,14,l)= 32.d0*c1*c2*c3*(four*c2-one)

              shl(1,15,l)= 32.d0*c2*(four*c3-one)*c3
     &                   -32.d0*c1*c2*(four*c3-one)-128.d0*c1*c2*c3
              shl(2,15,l)= 32.d0*c1*(four*c3-one)*c3
     &                   -32.d0*c1*c2*(four*c3-one)-128.d0*c1*c2*c3
              shl(3,15,l)= 32.d0*c1*c2*c3*(four*c3-one)

            end if
!           interpolação quíntica (p=5 e nen=21)(2012-10-24)
            if(nen.eq.21) then
              
        shl(1,1,l)= 
     &     (five/24.d0)*(five*c1-three)*(five*c1-two)*(five*c1-one)*c1
     &    +(five/24.d0)*(five*c1-four)*(five*c1-two)*(five*c1-one)*c1
     &    +(five/24.d0)*(five*c1-four)*(five*c1-three)*(five*c1-one)*c1
     &    +(five/24.d0)*(five*c1-four)*(five*c1-three)*(five*c1-two)*c1
     &    +(one/24.d0)*(five*c1-four)*(five*c1-three)
     &                *(five*c1-two)*(five*c1-one)
        shl(2,1,l)= zero
        shl(3,1,l)= (one/24.d0)*(five*c1-four)*(five*c1-three)
     &                         *(five*c1-two)*(five*c1-one)*c1
              
        shl(1,2,l)= zero
        shl(2,2,l)= 
     &     (five/24.d0)*(five*c2-three)*(five*c2-two)*(five*c2-one)*c2
     &    +(five/24.d0)*(five*c2-four)*(five*c2-two)*(five*c2-one)*c2
     &    +(five/24.d0)*(five*c2-four)*(five*c2-three)*(five*c2-one)*c2
     &    +(five/24.d0)*(five*c2-four)*(five*c2-three)*(five*c2-two)*c2
     &    +(one/24.d0)*(five*c2-four)*(five*c2-three)
     &                *(five*c2-two)*(five*c2-one)
        shl(3,2,l)= (one/24.d0)*(five*c2-four)*(five*c2-three)
     &                         *(five*c2-two)*(five*c2-one)*c2
              
        shl(1,3,l)= 
     &    -(five/24.d0)*(five*c3-three)*(five*c3-two)*(five*c3-one)*c3
     &    -(five/24.d0)*(five*c3-four)*(five*c3-two)*(five*c3-one)*c3
     &    -(five/24.d0)*(five*c3-four)*(five*c3-three)*(five*c3-one)*c3
     &    -(five/24.d0)*(five*c3-four)*(five*c3-three)*(five*c3-two)*c3
     &    -(one/24.d0)*(five*c3-four)*(five*c3-three)
     &                *(five*c3-two)*(five*c3-one)
        shl(2,3,l)= 
     &    -(five/24.d0)*(five*c3-three)*(five*c3-two)*(five*c3-one)*c3
     &    -(five/24.d0)*(five*c3-four)*(five*c3-two)*(five*c3-one)*c3
     &    -(five/24.d0)*(five*c3-four)*(five*c3-three)*(five*c3-one)*c3
     &    -(five/24.d0)*(five*c3-four)*(five*c3-three)*(five*c3-two)*c3
     &    -(one/24.d0)*(five*c3-four)*(five*c3-three)
     &                *(five*c3-two)*(five*c3-one)
        shl(3,3,l)= (one/24.d0)*(five*c3-four)*(five*c3-three)
     &                         *(five*c3-two)*(five*c3-one)*c3
              
        shl(1,4,l)= (125.d0/24.d0)*(five*c1-two)*(five*c1-one)*c1*c2
     &   + (125.d0/24.d0)*(five*c1-three)*(five*c1-one)*c1*c2
     &   + (125.d0/24.d0)*(five*c1-three)*(five*c1-two)*c1*c2
     &   + (25.d0/24.d0)*(five*c1-three)*(five*c1-two)
     &                  *(five*c1-one)*c2
        shl(2,4,l)= (25.d0/24.d0)*(five*c1-three)
     &                           *(five*c1-two)*(five*c1-one)*c1
        shl(3,4,l)= (25.d0/24.d0)*(five*c1-three)
     &                           *(five*c1-two)*(five*c1-one)*c1*c2
              
        shl(1,5,l)= (125.d0/12.d0)*(five*c1-one)*c1*(five*c2-one)*c2
     &    +(125.d0/12.d0)*(five*c1-two)*c1*(five*c2-one)*c2
     &    +(25.d0/12.d0)*(five*c1-two)*(five*c1-one)*(five*c2-one)*c2
        shl(2,5,l)= (125.d0/12.d0)*(five*c1-two)*(five*c1-one)*c1*c2
     &    +(25.d0/12.d0)*(five*c1-two)*(five*c1-one)*c1*(five*c2-one)
        shl(3,5,l)= (25.d0/12.d0)*(five*c1-two)
     &                           *(five*c1-one)*c1*(five*c2-one)*c2
              
        shl(1,6,l)= (125.d0/12.d0)*c1*(five*c2-two)*(five*c2-one)*c2
     &   + (25.d0/12.d0)*(five*c1-one)*(five*c2-two)*(five*c2-one)*c2
        shl(2,6,l)= (125.d0/12.d0)*(five*c1-one)*c1*(five*c2-one)*c2
     &   + (125.d0/12.d0)*(five*c1-one)*c1*(five*c2-two)*c2
     &   + (25.d0/12.d0)*(five*c1-one)*c1*(five*c2-two)*(five*c2-one)
        shl(3,6,l)= (25.d0/12.d0)*(five*c1-one)
     &                           *c1*(five*c2-two)*(five*c2-one)*c2

        shl(1,7,l)= (25.d0/24.d0)*(five*c2-three)
     &                           *(five*c2-two)*(five*c2-one)*c2
        shl(2,7,l)= (125.d0/24.d0)*c1*(five*c2-two)*(five*c2-one)*c2
     &    +(125.d0/24.d0)*c1*(five*c2-three)*(five*c2-one)*c2
     &    +(125.d0/24.d0)*c1*(five*c2-three)*(five*c2-two)*c2
     &    +(25.d0/24.d0)*c1*(five*c2-three)*(five*c2-two)*(five*c2-one)
        shl(3,7,l)= (25.d0/24.d0)*c1*
     &               (five*c2-three)*(five*c2-two)*(five*c2-one)*c2

        shl(1,8,l)= -(25.d0/24.d0)*(five*c2-three)
     &                            *(five*c2-two)*(five*c2-one)*c2
        shl(2,8,l)= 
     &    -(25.d0/24.d0)*(five*c2-three)*(five*c2-two)*(five*c2-one)*c2
     &    +(125.d0/24.d0)*c3*(five*c2-two)*(five*c2-one)*c2
     &    +(125.d0/24.d0)*c3*(five*c2-three)*(five*c2-one)*c2
     &    +(125.d0/24.d0)*c3*(five*c2-three)*(five*c2-two)*c2
     &    +(25.d0/24.d0)*c3*(five*c2-three)*(five*c2-two)*(five*c2-one)
        shl(3,8,l)= (25.d0/24.d0)*c3*
     &                 (five*c2-three)*(five*c2-two)*(five*c2-one)*c2

        shl(1,9,l)= - (125.d0/12.d0)*c3*(five*c2-two)*(five*c2-one)*c2
     &   - (25.d0/12.d0)*(five*c3-one)*(five*c2-two)*(five*c2-one)*c2
        shl(2,9,l)= -((125.d0/12.d0)*c3)*(five*c2-two)*(five*c2-one)*c2
     &    -((25.d0/12.d0)*(five*c3-one))*(five*c2-two)*(five*c2-one)*c2
     &    +((125.d0/12.d0)*(five*c3-one))*c3*(five*c2-one)*c2
     &    +((125.d0/12.d0)*(five*c3-one))*c3*(five*c2-two)*c2
     &    +((25.d0/12.d0)*(five*c3-one))*c3*(five*c2-two)*(five*c2-one)
        shl(3,9,l)= ((25.d0/12.d0)*(five*c3-one))
     &                     *c3*(five*c2-two)*(five*c2-one)*c2

        shl(1,10,l)= -((125.d0/12.d0)*(five*c3-one))*c3*(five*c2-one)*c2
     &    -((125.d0/12.d0)*(five*c3-two))*c3*(five*c2-one)*c2
     &    -((25.d0/12.d0)*(five*c3-two))*(five*c3-one)*(five*c2-one)*c2
        shl(2,10,l)= -((125.d0/12.d0)*(five*c3-one))*c3*(five*c2-one)*c2
     &    -((125.d0/12.d0)*(five*c3-two))*c3*(five*c2-one)*c2
     &    -((25.d0/12.d0)*(five*c3-two))*(five*c3-one)*(five*c2-one)*c2
     &    +((125.d0/12.d0)*(five*c3-two))*(five*c3-one)*c3*c2
     &    +((25.d0/12.d0)*(five*c3-two))*(five*c3-one)*c3*(five*c2-one)
        shl(3,10,l)= ((25.d0/12.d0)*(five*c3-two))*(five*c3-one)
     &                      *c3*(five*c2-one)*c2

        shl(1,11,l)= -((125.d0/24.d0)*(five*c3-two))*(five*c3-one)*c3*c2
     &   -((125.d0/24.d0)*(five*c3-three))*(five*c3-one)*c3*c2
     &   -((125.d0/24.d0)*(five*c3-three))*(five*c3-two)*c3*c2
     &   -((25.d0/24.d0)*(five*c3-three))*(five*c3-two)*(five*c3-one)*c2
        shl(2,11,l)= -((125.d0/24.d0)*(five*c3-two))*(five*c3-one)*c3*c2
     &    -((125.d0/24.d0)*(five*c3-three))*(five*c3-one)*c3*c2
     &    -((125.d0/24.d0)*(five*c3-three))*(five*c3-two)*c3*c2
     &  -((25.d0/24.d0)*(five*c3-three))*(five*c3-two)*(five*c3-one)*c2
     &  +((25.d0/24.d0)*(five*c3-three))*(five*c3-two)*(five*c3-one)*c3
        shl(3,11,l)= ((25.d0/24.d0)*(five*c3-three))
     &                     *(five*c3-two)*(five*c3-one)*c3*c2
              
        shl(1,12,l)= -((125.d0/24.d0)*(five*c3-two))*(five*c3-one)*c3*c1
     &    -((125.d0/24.d0)*(five*c3-three))*(five*c3-one)*c3*c1
     &    -((125.d0/24.d0)*(five*c3-three))*(five*c3-two)*c3*c1
     &  -((25.d0/24.d0)*(five*c3-three))*(five*c3-two)*(five*c3-one)*c1
     &  +((25.d0/24.d0)*(five*c3-three))*(five*c3-two)*(five*c3-one)*c3
        shl(2,12,l)= -((125.d0/24.d0)*(five*c3-two))*(five*c3-one)*c3*c1
     &    -((125.d0/24.d0)*(five*c3-three))*(five*c3-one)*c3*c1
     &    -((125.d0/24.d0)*(five*c3-three))*(five*c3-two)*c3*c1
     &   -((25.d0/24.d0)*(five*c3-three))*(five*c3-two)*(five*c3-one)*c1
        shl(3,12,l)= ((25.d0/24.d0)*(five*c3-three))
     &                   *(five*c3-two)*(five*c3-one)*c3*c1
              
        shl(1,13,l)= -((125.d0/12.d0)*(five*c3-one))*c3*(five*c1-one)*c1
     &    -((125.d0/12.d0)*(five*c3 - two))*c3*(five*c1 - one)*c1
     &    -((25.d0/12.d0)*(five*c3-two))*(five*c3-one)*(five*c1-one)*c1
     &    +((125.d0/12.d0)*(five*c3 - two))*(five*c3 - one)*c3*c1
     &    +((25.d0/12.d0)*(five*c3-two))*(five*c3-one)*c3*(five*c1-one)
        shl(2,13,l)= -((125.d0/12.d0)*(five*c3-one))*c3*(five*c1-one)*c1
     &    -((125.d0/12.d0)*(five*c3-two))*c3*(five*c1-one)*c1
     &    -((25.d0/12.d0)*(five*c3-two))*(five*c3-one)*(five*c1-one)*c1
        shl(3,13,l)= ((25.d0/12.d0)*(five*c3-two))
     &                      *(five*c3-one)*c3*(five*c1-one)*c1
              
        shl(1,14,l)= -((125.d0/12.d0)*c3)*(five*c1-two)*(five*c1-one)*c1
     &    -((25.d0/12.d0)*(five*c3-one))*(five*c1-two)*(five*c1-one)*c1
     &    +((125.d0/12.d0)*(five*c3-one))*c3*(five*c1-one)*c1
     &    +((125.d0/12.d0)*(five*c3-one))*c3*(five*c1-two)*c1
     &    +((25.d0/12.d0)*(five*c3-one))*c3*(five*c1-two)*(five*c1-one)
        shl(2,14,l)= -((125.d0/12.d0)*c3)*(five*c1-two)*(five*c1-one)*c1
     &    -((25.d0/12.d0)*(five*c3-one))*(five*c1-two)*(five*c1-one)*c1
        shl(3,14,l)= ((25.d0/12.d0)*(five*c3-one))
     &                       *c3*(five*c1-two)*(five*c1-one)*c1

        shl(1,15,l)= 
     &   -((25.d0/24.d0)*(five*c1-three))*(five*c1-two)*(five*c1-one)*c1
     &   +((125.d0/24.d0)*c3)*(five*c1-two)*(five*c1-one)*c1
     &   +((125.d0/24.d0)*c3)*(five*c1 - three)*(five*c1 - one)*c1
     &   +((125.d0/24.d0)*c3)*(five*c1 - three)*(five*c1 - two)*c1
     &   +((25.d0/24.d0)*c3)*(five*c1-three)*(five*c1-two)*(five*c1-one)
        shl(2,15,l)= -((25.d0/24.d0)*(five*c1 - three))
     &                       *(five*c1 - two)*(five*c1 - one)*c1
        shl(3,15,l)= (25.d0/24.d0)*c3*(five*c1-three)*(five*c1-two)
     &         *(five*c1-one)*c1

              
        shl(1,16,l)= ((625.d0/6.d0)*(five*c1-one))*c1*c2*c3
     &    +((625.d0/6.d0)*(five*c1 - two))*c1*c2*c3
     &    +((125.d0/6.d0)*(five*c1 - two))*(five*c1 - one)*c2*c3
     &    -((125.d0/6.d0)*(five*c1 - two))*(five*c1 - one)*c1*c2
        shl(2,16,l)= ((125.d0/6.d0)*c3)*(five*c1-two)*(five*c1-one)*c1
     &    -((125.d0/6.d0)*(five*c1 - two))*(five*c1-one)*c1*c2
        shl(3,16,l)=((125.d0/6.d0)*(five*c1-two))*(five*c1-one)*c1*c2*c3

        shl(1,17,l)= ((625.d0/4.d0)*(five*c2 - one))*c1*c2*c3
     &    +((125.d0/4.d0)*(five*c1 - one))*(five*c2 - one)*c2*c3
     &    -((125.d0/4.d0)*(five*c1 - one))*c1*(five*c2 - one)*c2
        shl(2,17,l)= ((625.d0/4.d0)*(five*c1 - one))*c1*c2*c3
     &    +((125.d0/4.d0)*(five*c1 - one))*(five*c2 - one)*c1*c3
     &    -((125.d0/4.d0)*(five*c1 - one))*c1*(five*c2 - one)*c2
       shl(3,17,l)=((125.d0/4.d0)*(five*c1-one))*(five*c2-one)*c1*c2*c3
              
        shl(1,18,l)= ((125.d0/6.d0)*c3)*(five*c2-two)*(five*c2 - one)*c2
     &    -((125.d0/6.d0))*c1*(five*c2 - two)*(five*c2 - one)*c2
        shl(2,18,l)= ((625.d0/6.d0)*(five*c2 - one))*c1*c2*c3
     &    +((625.d0/6.d0)*(five*c2 - two))*c1*c2*c3
     &    +((125.d0/6.d0)*(five*c2 - two))*(five*c2 - one)*c1*c3
     &    -((125.d0/6.d0))*c1*(five*c2 - two)*(five*c2 - one)*c2
        shl(3,18,l)=((125.d0/6.d0)*(five*c2-two))*(five*c2-one)*c1*c2*c3
              
        shl(1,19,l)= -((625.d0/4.d0)*(five*c2 - one))*c1*c2*c3
     &    +((125.d0/4.d0)*(five*c3 - one))*c3*(five*c2 - one)*c2
     &    -((125.d0/4.d0)*(five*c2 - one))*(five*c3 - one)*c1*c2
        shl(2,19,l)= ((625.d0/4.d0)*(five*c3 - one))*c1*c2*c3
     &    -((625.d0/4.d0)*(five*c2 - one))*c1*c2*c3
     &    +((125.d0/4.d0)*(five*c2 - one))*(five*c3 - one)*c1*c3
     &    -((125.d0/4.d0)*(five*c2 - one))*(five*c3 - one)*c1*c2
       shl(3,19,l)= ((125.d0/4.d0)*(five*c2-one))*(five*c3-one)*c1*c2*c3
              
        shl(1,20,l)= -((625.d0/6.d0)*(five*c3 - one))*c1*c2*c3
     &    -((625.d0/6.d0)*(five*c3 - two))*c1*c2*c3
     &    +((125.d0/6.d0)*(five*c3 - two))*(five*c3 - one)*c3*c2
     &    -((125.d0/6.d0)*(five*c3 - two))*(five*c3 - one)*c1*c2
        shl(2,20,l)= -((625.d0/6.d0)*(five*c3 - one))*c1*c2*c3
     &    -((625.d0/6.d0)*(five*c3 - two))*c1*c2*c3
     &    +((125.d0/6.d0)*(five*c3 - two))*(five*c3 - one)*c3*c1
     &    -((125.d0/6.d0)*(five*c3 - two))*(five*c3 - one)*c1*c2
       shl(3,20,l)= ((125.d0/6.d0)*(five*c3-two))*(five*c3-one)*c1*c2*c3

        shl(1,21,l)= ((625.d0/4.d0)*(five*c3 - one))*c1*c2*c3
     &    -((625.d0/4.d0)*(five*c1 - one))*c1*c2*c3
     &    +((125.d0/4.d0)*(five*c1 - one))*(five*c3 - one)*c2*c3
     &    -((125.d0/4.d0)*(five*c1 - one))*(five*c3 - one)*c1*c2
        shl(2,21,l)= -((625.d0/4.d0)*(five*c1 - one))*c1*c2*c3
     &    +((125.d0/4.d0)*(five*c3 - one))*c3*(five*c1 - one)*c1
     &    -((125.d0/4.d0)*(five*c1 - one))*(five*c3 - one)*c1*c2
       shl(3,21,l)= ((125.d0/4.d0)*(five*c1-one))*(five*c3-one)*c1*c2*c3

            end if
            
  200      continue 
c 
           return 
           end
C 
c******* nova subrotina ***************************************************************
c                                              
      subroutine shltpbk(shl,nen,nside,nints)
c
c     Objetivo: calcular pesos, funcoes de interpolacao e derivadas locais
c               para elementos triangulares
c
c----------------------------------------------------------------------
c
c 
c.... program to calculate integration-rule weights, shape functions 
c        and local derivatives for a triangular element 
c 
c        c1, c2, c3 = local element coordinates ("l1", "l2", "l3".) 
c        shl(j,i,l) = local ("j") derivative of shape function 
c        shl(3,i,l) = local  shape function 
c              w(l) = integration-rule weight 
c                 i = local node number 
c                 l = integration point number 
c              nint = number of integration points, eq. 1 or 4 
c 
      implicit real*8 (a-h,o-z) 
c 
c.... remove above card for single precision operation 
c 
      dimension shl(3,nen,*),cl1(18),cl2(18),cl3(18),csi(6),ra(6)
      data   zero,pt1667,pt25,pt5 
     &      /0.0d0,0.1666666666666667d0,0.25d0,0.5d0/, 
     &       one,two,three,four,five,six 
     &      /1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0/ 
c 
      pt3s2 = three/two 
      pt9s2 = 9.d0/two
      pt27s2= 27.d0/two
      pt2s3 = two/three
      pt1s6 = one/six
      pt32s3= 32.d0/three
      pt8s3 = 8.d0/three
c
      if (nints.eq.1) then
         ra(1) = zero
      endif
c
      if (nints.eq.2) then
         ra(1)=-.577350269189626
         ra(2)= .577350269189625
      endif
c
      if (nints.eq.3) then
         ra(1)=-.774596669241483
         ra(2)= zero
         ra(3)= .774596669241483
      endif
c
      if (nints.eq.4) then
         ra(1)=-.861136311594053
         ra(2)=-.339981043584856
         ra(3)= .339981043584856
         ra(4)= .861136311594053
      endif
c
       if(nints.eq.5) then
        ra(1)=-.906179845938664
        ra(2)=-.538469310105683
        ra(3)= zero
        ra(4)= .538469310105683
        ra(5)= .906179845938664
       endif
c
       if(nints.eq.6) then
         ra(1)=-.932469514203152
         ra(2)=-.661209386466265
         ra(3)=-.238619186083197
         ra(4)= .238619186083197
         ra(5)= .661209386466365
         ra(6)= .932469514203152
        endif
c
       if(nints.eq.7) then
         ra(1)=-.949107912342759
         ra(2)=-.741531185599394
         ra(3)=-.405845151377397
         ra(4)= zero
         ra(5)= .405845151377397
         ra(6)= .741531185599394
         ra(7)= .949107912342759
        endif

c
      lb = 0
      do 100 ns=1,3
      do 200 ls=1,nints
c
      lb = lb + 1
c
      if(ns.eq.1) then
        cl3(ls) = 0.d00
        cl1(ls) = (1.d00 - ra(ls))/two
        cl2(ls) = (1.d00 + ra(ls))/two
      end if 
c
      if(ns.eq.2) then
         cl1(ls) = 0.d00
         cl2(ls) = (1.d00 - ra(ls))/two
         cl3(ls) = (1.d00 + ra(ls))/two
      end if 
c
      if(ns.eq.3) then
         cl2(ls) = 0.d00
         cl3(ls) = (1.d00 - ra(ls))/two
         cl1(ls) = (1.d00 + ra(ls))/two
      end if 
c 
            c1 = cl1(ls) 
            c2 = cl2(ls) 
            c3 = cl3(ls) 
c 
            if(nen.eq.1) then
              shl(1,1,lb)= zero
              shl(2,1,lb)= zero
              shl(3,1,lb)= one
            end if             
!          acrescentei o if (2012-10-08)
!           interpolação linear (p=1 e nen=3)(2012-10-08)
            if(nen.eq.3) then
              shl(1,1,lb)= one
              shl(2,1,lb)= zero
              shl(3,1,lb)= c1
              shl(1,2,lb)= zero
              shl(2,2,lb)= one
              shl(3,2,lb)= c2
              shl(1,3,lb)=-one
              shl(2,3,lb)=-one
              shl(3,3,lb)= c3
            end if 
!           interpolação quadrática (p=2 e nen=6)(2012-10-08)
            if(nen.eq.6) then
              shl(1,1,lb)= four*c1-one
              shl(2,1,lb)= zero
              shl(3,1,lb)= (two*c1 - one)*c1
              shl(1,2,lb)= zero
              shl(2,2,lb)= four*c2-one
              shl(3,2,lb)= (two*c2 - one)*c2
              shl(1,3,lb)= one - four*c3
              shl(2,3,lb)= one - four*c3
              shl(3,3,lb)= (two*c3 - one)*c3
              shl(1,4,lb)= four * c2
              shl(2,4,lb)= four * c1
              shl(3,4,lb)= four * c1 * c2
              shl(1,5,lb)=-four * c2
              shl(2,5,lb)= four * (c3 - c2)
              shl(3,5,lb)= four * c2 * c3
              shl(1,6,lb)= four * (c3 - c1)
              shl(2,6,lb)=-four * c1
              shl(3,6,lb)= four * c3 * c1
            end if
c
            if (nen.eq.10) then 
 
                  shl(1,1,lb)= pt5*((three*c1-one)*(three*c1-two) 
     &            +c1*three*(three*c1-two)+c1*three*(three*c1-one)) 
                  shl(2,1,lb)= zero 
                  shl(3,1,lb)= pt5*c1*(three*c1-one)*(three*c1-two) 
 
                  shl(1,2,lb)=zero 
                  shl(2,2,lb)=  pt5*((three*c2-one)*(three*c2-two) 
     &            +c2*three*(three*c2-two)+c2*three*(three*c2-one)) 
                  shl(3,2,lb)=  pt5*c2*(three*c2-one)*(three*c2-two) 
 
                  shl(1,3,lb)= -pt5*((three*c3-one)*(three*c3-two) 
     &            +c3*three*(three*c3-two)+c3*three*(three*c3-one)) 
                  shl(2,3,lb)=-pt5*((three*c3-one)*(three*c3-two) 
     &            +c3*three*(three*c3-two)+c3*three*(three*c3-one)) 
                  shl(3,3,lb)=  pt5*c3*(three*c3-one)*(three*c3-two) 
 
 
                  shl(1,4,lb) = (9.0d0/two)*(c2*(three*c1-one) 
     &            +c1*c2*three) 
                  shl(2,4,lb) = (9.0d0/two)*c1*(three*c1-one) 
                  shl(3,4,lb) = (9.0d0/two)*c1*c2*(three*c1-one) 
 
                  shl(1,5,lb) = (9.0d0/two)*c2*(three*c2-one) 
                  shl(2,5,lb) = (9.0d0/two)*(c1*(three*c2-one) 
     &            +c1*c2*three) 
                  shl(3,5,lb) = (9.0d0/two)*c1*c2*(three*c2-one) 
 
 
                 shl(1,6,lb) =- (9.0d0/two)*c2*(three*c2-one) 
                 shl(2,6,lb) =  (9.0d0/two)*(c3*(three*c2-one) 
     &           +c2*c3*three-c2*(three*c2-one)) 
                 shl(3,6,lb) = (9.0d0/two)*c3*c2*(three*c2-one) 
 
 
                  shl(1,7,lb) = (9.0d0/two)*(-c2*(three*c3-one) 
     &            -c2*c3*three) 
                  shl(2,7,lb) =  (9.0d0/two)*(c3*(three*c3-one) 
     &            -c2*c3*three -c2*(three*c3-one)) 
                  shl(3,7,lb) = (9.0d0/two)*c2*c3*(three*c3-one) 
 
                  shl(1,9,lb)= (9.0d0/two)*(-c1*(three*c1-one) 
     &            +c1*c3*three +c3*(three*c1-one)) 
                  shl(2,9,lb) =  (9.0d0/two)*(-c1*(three*c1-one)) 
                  shl(3,9,lb) = (9.0d0/two)*c3*c1*(three*c1-one) 
 
 
                  shl(1,8,lb) = (9.0d0/two)*(-c1*(three*c3-one) 
     &            -c1*c3*three +c3*(three*c3-one)) 
                  shl(2,8,lb) =  (9.0d0/two)*( 
     &            -c1*c3*three -c1*(three*c3-one)) 
                  shl(3,8,lb) = (9.0d0/two)*c1*c3*(three*c3-one) 
 
                 shl(1,10,lb) = 27.0d0*(c2*c3-c1*c2) 
                 shl(2,10,lb) = 27.0d0*(c1*c3-c1*c2) 
                 shl(3,10,lb) = 27.0d0*c1*c2*c3 
          endif 
!           interpolação quártica (p=4 e nen=15)(2012-10-23)
            if(nen.eq.15) then
           
              shl(1,1,lb)= pt2s3*(four*c1-two)*(four*c1-one)*c1
     &                  + pt2s3*(four*c1-three)*(four*c1-one)*c1
     &                  + pt2s3*(four*c1-three)*(four*c1-two)*c1
     &            + pt1s6*(four*c1-three)*(four*c1-two)*(four*c1-one)
              shl(2,1,lb)= zero
              shl(3,1,lb)= pt1s6*(four*c1-three)*(four*c1-two)*
     &                  (four*c1-one)*c1
              
              shl(1,2,lb)= zero
              shl(2,2,lb)= pt2s3*(four*c2-two)*(four*c2-one)*c2
     &                  + pt2s3*(four*c2-three)*(four*c2-one)*c2
     &                  + pt2s3*(four*c2-three)*(four*c2-two)*c2
     &            + pt1s6*(four*c2-three)*(four*c2-two)*(four*c2-one)
              shl(3,2,lb)= pt1s6*(four*c2-three)*(four*c2-two)*
     &                  (four*c2-one)*c2
              
              shl(1,3,lb)= -pt2s3*(four*c3-two)*(four*c3-one)*c3
     &                  -pt2s3*(four*c3-three)*(four*c3-one)*c3
     &                  -pt2s3*(four*c3-three)*(four*c3-two)*c3
     &             -pt1s6*(four*c3-three)*(four*c3-two)*(four*c3-one)
              shl(2,3,lb)= -pt2s3*(four*c3-two)*(four*c3-one)*c3
     &                  -pt2s3*(four*c3-three)*(four*c3-one)*c3
     &                  -pt2s3*(four*c3-three)*(four*c3-two)*c3
     &             -pt1s6*(four*c3-three)*(four*c3-two)*(four*c3-one)
              shl(3,3,lb)= pt1s6*(four*c3-three)*(four*c3-two)*
     &                   (four*c3-one)*c3
              
              shl(1,4,lb)= pt32s3*c2*(four*c1-one)*c1
     &                  +pt32s3*c2*(four*c1-two)*c1
     &                  +pt8s3*c2*(four*c1-two)*(four*c1-one)
              shl(2,4,lb)= pt8s3*(four*c1-two)*(four*c1-one)*c1
              shl(3,4,lb)= pt8s3*c2*(four*c1-two)*(four*c1-one)*c1
              
              shl(1,5,lb)= 16.d0*c1*(four*c2-one)*c2 
     &                  + four*(four*c1-one)*(four*c2-one)*c2
              shl(2,5,lb)= 16.d0*c2*(four*c1-one)*c1
     &                  + four*(four*c1-one)*c1*(four*c2-one)
              shl(3,5,lb)= four*(four*c1-one)*c1*(four*c2-one)*c2
              
              shl(1,6,lb)= pt8s3*(four*c2-two)*(four*c2-one)*c2
              shl(2,6,lb)= pt32s3*c1*(four*c2-one)*c2 
     &                  + pt32s3*c1*(four*c2-two)*c2 
     &                  + pt8s3*c1*(four*c2-two)*(four*c2-one)
              shl(3,6,lb)= pt8s3*c1*(four*c2-two)*(four*c2-one)*c2

              shl(1,7,lb)= -pt8s3*(four*c2-two)*(four*c2-one)*c2
              shl(2,7,lb)= -pt8s3*(four*c2-two)*(four*c2-one)*c2 
     &                  + pt32s3*c3*(four*c2-one)*c2
     &                  + pt32s3*c3*(four*c2-two)*c2
     &                  + pt8s3*c3*(four*c2-two)*(four*c2-one)
              shl(3,7,lb)= pt8s3*c3*(four*c2-two)*(four*c2-one)*c2

              shl(1,8,lb)= -16.d0*(four*c2-one)*c2*c3
     &                  -four*(four*c2-one)*c2*(four*c3-one)
              shl(2,8,lb)= 16.d0*c2*(four*c3-one)*c3
     &                  +four*(four*c2-one)*(four*c3-one)*c3
     &                  -16.d0*(four*c2-one)*c2*c3
     &                  -four*(four*c2-one)*c2*(four*c3-one)
              shl(3,8,lb)= four*(four*c2-one)*c2*(four*c3-one)*c3

              shl(1,9,lb)= -pt32s3*c2*(four*c3-one)*c3
     &                  -pt32s3*c2*(four*c3-two)*c3
     &                  -pt8s3*c2*(four*c3-two)*(four*c3-one)
              shl(2,9,lb)= pt8s3*(four*c3-two)*(four*c3-one)*c3
     &                  -pt32s3*c2*(four*c3-one)*c3
     &                  -pt32s3*c2*(four*c3-two)*c3
     &                  -pt8s3*c2*(four*c3-two)*(four*c3-one)
              shl(3,9,lb)= pt8s3*c2*(four*c3-two)*(four*c3-one)*c3

              shl(1,10,lb)= pt8s3*(four*c3-two)*(four*c3-one)*c3
     &                   -pt32s3*c1*(four*c3-one)*c3
     &                   -pt32s3*c1*(four*c3-two)*c3
     &                   -pt8s3*c1*(four*c3-two)*(four*c3-one)
              shl(2,10,lb)= -pt32s3*c1*(four*c3-one)*c3
     &                   -pt32s3*c1*(four*c3-two)*c3
     &                   -pt8s3*c1*(four*c3-two)*(four*c3-one)
              shl(3,10,lb)= pt8s3*c1*(four*c3-two)*(four*c3-one)*c3

              shl(1,11,lb)= -16.d0*c3*(four*c1-one)*c1
     &                   -4.d0*(four*c3-one)*(four*c1-one)*c1
     &                   +16.d0*c1*(four*c3-one)*c3
     &                   +4.d0*(four*c3-one)*c3*(four*c1-one)
              shl(2,11,lb)= -16.d0*c3*(four*c1-one)*c1
     &                   -4.d0*(four*c3-one)*(four*c1-one)*c1
              shl(3,11,lb)= 4.d0*(four*c3-one)*c3*(four*c1-one)*c1
              
              shl(1,12,lb)= -pt8s3*(four*c1-two)*(four*c1-one)*c1
     &                   +pt32s3*c3*(four*c1-one)*c1
     &                   +pt32s3*c3*(four*c1-two)*c1
     &                   +pt8s3*c3*(four*c1-two)*(four*c1-one)
              shl(2,12,lb)= -pt8s3*(four*c1-two)*(four*c1-one)*c1
              shl(3,12,lb)= pt8s3*c3*(four*c1-two)*(four*c1-one)*c1
              
              shl(1,13,lb)= 32.d0*c2*c3*(four*c1-one)
     &                   -32.d0*c1*c2*(four*c1-one)+128.d0*c1*c2*c3
              shl(2,13,lb)= 32.d0*c1*c3*(four*c1-one)
     &                   -32.d0*c1*c2*(four*c1-one)
              shl(3,13,lb)= 32.d0*c1*c2*c3*(four*c1-one)
              
              shl(1,14,lb)= 32.d0*(four*c2-one)*c2*c3
     &                   -32.d0*c1*(four*c2-one)*c2
              shl(2,14,lb)= 32.d0*c1*c3*(four*c2-one)
     &                   -32.d0*c1*c2*(four*c2-one)+128.d0*c1*c2*c3
              shl(3,14,lb)= 32.d0*c1*c2*c3*(four*c2-one)

              shl(1,15,lb)= 32.d0*c2*(four*c3-one)*c3
     &                   -32.d0*c1*c2*(four*c3-one)-128.d0*c1*c2*c3
              shl(2,15,lb)= 32.d0*c1*(four*c3-one)*c3
     &                   -32.d0*c1*c2*(four*c3-one)-128.d0*c1*c2*c3
              shl(3,15,lb)= 32.d0*c1*c2*c3*(four*c3-one)

            end if
!           interpolação quíntica (p=5 e nen=21)(2012-10-24)
            if(nen.eq.21) then
              
        shl(1,1,lb)= 
     &     (five/24.d0)*(five*c1-three)*(five*c1-two)*(five*c1-one)*c1
     &    +(five/24.d0)*(five*c1-four)*(five*c1-two)*(five*c1-one)*c1
     &    +(five/24.d0)*(five*c1-four)*(five*c1-three)*(five*c1-one)*c1
     &    +(five/24.d0)*(five*c1-four)*(five*c1-three)*(five*c1-two)*c1
     &    +(one/24.d0)*(five*c1-four)*(five*c1-three)
     &                *(five*c1-two)*(five*c1-one)
        shl(2,1,lb)= zero
        shl(3,1,lb)= (one/24.d0)*(five*c1-four)*(five*c1-three)
     &                    *(five*c1-two)*(five*c1-one)*c1
              
        shl(1,2,lb)= zero
        shl(2,2,lb)= 
     &     (five/24.d0)*(five*c2-three)*(five*c2-two)*(five*c2-one)*c2
     &    +(five/24.d0)*(five*c2-four)*(five*c2-two)*(five*c2-one)*c2
     &    +(five/24.d0)*(five*c2-four)*(five*c2-three)*(five*c2-one)*c2
     &    +(five/24.d0)*(five*c2-four)*(five*c2-three)*(five*c2-two)*c2
     &    +(one/24.d0)*(five*c2-four)*(five*c2-three)
     &                *(five*c2-two)*(five*c2-one)
        shl(3,2,lb)= (one/24.d0)*(five*c2-four)*(five*c2-three)
     &                    *(five*c2-two)*(five*c2-one)*c2
              
        shl(1,3,lb)= 
     &    -(five/24.d0)*(five*c3-three)*(five*c3-two)*(five*c3-one)*c3
     &    -(five/24.d0)*(five*c3-four)*(five*c3-two)*(five*c3-one)*c3
     &    -(five/24.d0)*(five*c3-four)*(five*c3-three)*(five*c3-one)*c3
     &    -(five/24.d0)*(five*c3-four)*(five*c3-three)*(five*c3-two)*c3
     &    -(one/24.d0)*(five*c3-four)*(five*c3-three)
     &                 *(five*c3-two)*(five*c3-one)
        shl(2,3,lb)= 
     &    -(five/24.d0)*(five*c3-three)*(five*c3-two)*(five*c3-one)*c3
     &    -(five/24.d0)*(five*c3-four)*(five*c3-two)*(five*c3-one)*c3
     &    -(five/24.d0)*(five*c3-four)*(five*c3-three)*(five*c3-one)*c3
     &    -(five/24.d0)*(five*c3-four)*(five*c3-three)*(five*c3-two)*c3
     &    -(one/24.d0)*(five*c3-four)*(five*c3-three)
     &                  *(five*c3-two)*(five*c3-one)
        shl(3,3,lb)= (one/24.d0)*(five*c3-four)*(five*c3-three)
     &                    *(five*c3-two)*(five*c3-one)*c3
              
        shl(1,4,lb)= (125.d0/24.d0)*(five*c1-two)*(five*c1-one)*c1*c2
     &   + (125.d0/24.d0)*(five*c1-three)*(five*c1-one)*c1*c2
     &   + (125.d0/24.d0)*(five*c1-three)*(five*c1-two)*c1*c2
     &   + (25.d0/24.d0)*(five*c1-three)*(five*c1-two)
     &                  *(five*c1-one)*c2
        shl(2,4,lb)= (25.d0/24.d0)*(five*c1-three)
     &                    *(five*c1-two)*(five*c1-one)*c1
        shl(3,4,lb)= (25.d0/24.d0)*(five*c1-three)
     &                    *(five*c1-two)*(five*c1-one)*c1*c2
              
        shl(1,5,lb)= (125.d0/12.d0)*(five*c1-one)*c1*(five*c2-one)*c2
     &    +(125.d0/12.d0)*(five*c1-two)*c1*(five*c2-one)*c2
     &    +(25.d0/12.d0)*(five*c1-two)*(five*c1-one)*(five*c2-one)*c2
        shl(2,5,lb)= (125.d0/12.d0)*(five*c1-two)*(five*c1-one)*c1*c2
     &    +(25.d0/12.d0)*(five*c1-two)*(five*c1-one)*c1*(five*c2-one)
        shl(3,5,lb)= (25.d0/12.d0)*(five*c1-two)
     &                    *(five*c1-one)*c1*(five*c2-one)*c2
              
        shl(1,6,lb)= (125.d0/12.d0)*c1*(five*c2-two)*(five*c2-one)*c2
     &   + (25.d0/12.d0)*(five*c1-one)*(five*c2-two)*(five*c2-one)*c2
        shl(2,6,lb)= (125.d0/12.d0)*(five*c1-one)*c1*(five*c2-one)*c2
     &   + (125.d0/12.d0)*(five*c1-one)*c1*(five*c2-two)*c2
     &   + (25.d0/12.d0)*(five*c1-one)*c1*(five*c2-two)*(five*c2-one)
        shl(3,6,lb)= (25.d0/12.d0)*(five*c1-one)
     &                    *c1*(five*c2-two)*(five*c2-one)*c2

        shl(1,7,lb)= (25.d0/24.d0)*(five*c2-three)
     &                    *(five*c2-two)*(five*c2-one)*c2
        shl(2,7,lb)= (125.d0/24.d0)*c1*(five*c2-two)*(five*c2-one)*c2
     &    +(125.d0/24.d0)*c1*(five*c2-three)*(five*c2-one)*c2
     &    +(125.d0/24.d0)*c1*(five*c2-three)*(five*c2-two)*c2
     &    +(25.d0/24.d0)*c1*(five*c2-three)*(five*c2-two)*(five*c2-one)
        shl(3,7,lb)= (25.d0/24.d0)*c1*
     &                (five*c2-three)*(five*c2-two)*(five*c2-one)*c2

        shl(1,8,lb)= -(25.d0/24.d0)*(five*c2-three)
     &                     *(five*c2-two)*(five*c2-one)*c2
        shl(2,8,lb)= 
     &    -(25.d0/24.d0)*(five*c2-three)*(five*c2-two)*(five*c2-one)*c2
     &    +(125.d0/24.d0)*c3*(five*c2-two)*(five*c2-one)*c2
     &    +(125.d0/24.d0)*c3*(five*c2-three)*(five*c2-one)*c2
     &    +(125.d0/24.d0)*c3*(five*c2-three)*(five*c2-two)*c2
     &    +(25.d0/24.d0)*c3*(five*c2-three)*(five*c2-two)*(five*c2-one)
        shl(3,8,lb)= (25.d0/24.d0)*c3*
     &                 (five*c2-three)*(five*c2-two)*(five*c2-one)*c2

        shl(1,9,lb)= - (125.d0/12.d0)*c3*(five*c2-two)*(five*c2-one)*c2
     &   - (25.d0/12.d0)*(five*c3-one)*(five*c2-two)*(five*c2-one)*c2
        shl(2,9,lb)= -((125.d0/12.d0)*c3)*(five*c2-two)*(five*c2-one)*c2
     &   -((25.d0/12.d0)*(five*c3-one))*(five*c2-two)*(five*c2-one)*c2
     &   +((125.d0/12.d0)*(five*c3-one))*c3*(five*c2-one)*c2
     &   +((125.d0/12.d0)*(five*c3-one))*c3*(five*c2-two)*c2
     &   +((25.d0/12.d0)*(five*c3-one))*c3*(five*c2-two)*(five*c2-one)
        shl(3,9,lb)= ((25.d0/12.d0)*(five*c3-one))
     &                     *c3*(five*c2-two)*(five*c2-one)*c2

       shl(1,10,lb)= -((125.d0/12.d0)*(five*c3-one))*c3*(five*c2-one)*c2
     &   -((125.d0/12.d0)*(five*c3-two))*c3*(five*c2-one)*c2
     &   -((25.d0/12.d0)*(five*c3-two))*(five*c3-one)*(five*c2-one)*c2
       shl(2,10,lb)= -((125.d0/12.d0)*(five*c3-one))*c3*(five*c2-one)*c2
     &    -((125.d0/12.d0)*(five*c3-two))*c3*(five*c2-one)*c2
     &    -((25.d0/12.d0)*(five*c3-two))*(five*c3-one)*(five*c2-one)*c2
     &    +((125.d0/12.d0)*(five*c3-two))*(five*c3-one)*c3*c2
     &    +((25.d0/12.d0)*(five*c3-two))*(five*c3-one)*c3*(five*c2-one)
        shl(3,10,lb)= ((25.d0/12.d0)*(five*c3-two))*(five*c3-one)
     &                      *c3*(five*c2-one)*c2

       shl(1,11,lb)= -((125.d0/24.d0)*(five*c3-two))*(five*c3-one)*c3*c2
     &   -((125.d0/24.d0)*(five*c3-three))*(five*c3-one)*c3*c2
     &   -((125.d0/24.d0)*(five*c3-three))*(five*c3-two)*c3*c2
     &  -((25.d0/24.d0)*(five*c3-three))*(five*c3-two)*(five*c3-one)*c2
       shl(2,11,lb)= -((125.d0/24.d0)*(five*c3-two))*(five*c3-one)*c3*c2
     &   -((125.d0/24.d0)*(five*c3-three))*(five*c3-one)*c3*c2
     &   -((125.d0/24.d0)*(five*c3-three))*(five*c3-two)*c3*c2
     &   -((25.d0/24.d0)*(five*c3-three))*(five*c3-two)*(five*c3-one)*c2
     &   +((25.d0/24.d0)*(five*c3-three))*(five*c3-two)*(five*c3-one)*c3
        shl(3,11,lb)= ((25.d0/24.d0)*(five*c3-three))
     &                     *(five*c3-two)*(five*c3-one)*c3*c2
              
       shl(1,12,lb)= -((125.d0/24.d0)*(five*c3-two))*(five*c3-one)*c3*c1
     &             -((125.d0/24.d0)*(five*c3-three))*(five*c3-one)*c3*c1
     &             -((125.d0/24.d0)*(five*c3-three))*(five*c3-two)*c3*c1
     &   -((25.d0/24.d0)*(five*c3-three))*(five*c3-two)*(five*c3-one)*c1
     &   +((25.d0/24.d0)*(five*c3-three))*(five*c3-two)*(five*c3-one)*c3
       shl(2,12,lb)= -((125.d0/24.d0)*(five*c3-two))*(five*c3-one)*c3*c1
     &             -((125.d0/24.d0)*(five*c3-three))*(five*c3-one)*c3*c1
     &             -((125.d0/24.d0)*(five*c3-three))*(five*c3-two)*c3*c1
     &   -((25.d0/24.d0)*(five*c3-three))*(five*c3-two)*(five*c3-one)*c1
              shl(3,12,lb)= ((25.d0/24.d0)*(five*c3-three))
     &                      *(five*c3-two)*(five*c3-one)*c3*c1
              
       shl(1,13,lb)= -((125.d0/12.d0)*(five*c3-one))*c3*(five*c1-one)*c1
     &           -((125.d0/12.d0)*(five*c3 - two))*c3*(five*c1 - one)*c1
     &     -((25.d0/12.d0)*(five*c3-two))*(five*c3-one)*(five*c1-one)*c1
     &           +((125.d0/12.d0)*(five*c3 - two))*(five*c3 - one)*c3*c1
     &     +((25.d0/12.d0)*(five*c3-two))*(five*c3-one)*c3*(five*c1-one)
       shl(2,13,lb)= -((125.d0/12.d0)*(five*c3-one))*c3*(five*c1-one)*c1
     &             -((125.d0/12.d0)*(five*c3-two))*c3*(five*c1-one)*c1
     &     -((25.d0/12.d0)*(five*c3-two))*(five*c3-one)*(five*c1-one)*c1
              shl(3,13,lb)= ((25.d0/12.d0)*(five*c3-two))
     &                      *(five*c3-one)*c3*(five*c1-one)*c1
              
       shl(1,14,lb)= -((125.d0/12.d0)*c3)*(five*c1-two)*(five*c1-one)*c1
     &     -((25.d0/12.d0)*(five*c3-one))*(five*c1-two)*(five*c1-one)*c1
     &             +((125.d0/12.d0)*(five*c3-one))*c3*(five*c1-one)*c1
     &             +((125.d0/12.d0)*(five*c3-one))*c3*(five*c1-two)*c1
     &     +((25.d0/12.d0)*(five*c3-one))*c3*(five*c1-two)*(five*c1-one)
       shl(2,14,lb)= -((125.d0/12.d0)*c3)*(five*c1-two)*(five*c1-one)*c1
     &     -((25.d0/12.d0)*(five*c3-one))*(five*c1-two)*(five*c1-one)*c1
              shl(3,14,lb)= ((25.d0/12.d0)*(five*c3-one))
     &                       *c3*(five*c1-two)*(five*c1-one)*c1

              shl(1,15,lb)= 
     &   -((25.d0/24.d0)*(five*c1-three))*(five*c1-two)*(five*c1-one)*c1
     &          +((125.d0/24.d0)*c3)*(five*c1-two)*(five*c1-one)*c1
     &        +((125.d0/24.d0)*c3)*(five*c1 - three)*(five*c1 - one)*c1
     &        +((125.d0/24.d0)*c3)*(five*c1 - three)*(five*c1 - two)*c1
     &   +((25.d0/24.d0)*c3)*(five*c1-three)*(five*c1-two)*(five*c1-one)

              shl(2,15,lb)= -((25.d0/24.d0)*(five*c1 - three))
     &                       *(five*c1 - two)*(five*c1 - one)*c1
            shl(3,15,lb)= (25.d0/24.d0)*c3*(five*c1-three)*(five*c1-two)
     &         *(five*c1-one)*c1

              
              shl(1,16,lb)= ((625.d0/6.d0)*(five*c1-one))*c1*c2*c3
     &                   +((625.d0/6.d0)*(five*c1 - two))*c1*c2*c3
     &           +((125.d0/6.d0)*(five*c1 - two))*(five*c1 - one)*c2*c3
     &           -((125.d0/6.d0)*(five*c1 - two))*(five*c1 - one)*c1*c2
        shl(2,16,lb)= ((125.d0/6.d0)*c3)*(five*c1-two)*(five*c1-one)*c1
     &              -((125.d0/6.d0)*(five*c1 - two))*(five*c1-one)*c1*c2
       shl(3,16,lb)=((125.d0/6.d0)*(five*c1-two))*(five*c1-one)*c1*c2*c3

              shl(1,17,lb)= ((625.d0/4.d0)*(five*c2 - one))*c1*c2*c3
     &           +((125.d0/4.d0)*(five*c1 - one))*(five*c2 - one)*c2*c3
     &            -((125.d0/4.d0)*(five*c1 - one))*c1*(five*c2 - one)*c2
              shl(2,17,lb)= ((625.d0/4.d0)*(five*c1 - one))*c1*c2*c3
     &           +((125.d0/4.d0)*(five*c1 - one))*(five*c2 - one)*c1*c3
     &           -((125.d0/4.d0)*(five*c1 - one))*c1*(five*c2 - one)*c2
       shl(3,17,lb)=((125.d0/4.d0)*(five*c1-one))*(five*c2-one)*c1*c2*c3
              
       shl(1,18,lb)= ((125.d0/6.d0)*c3)*(five*c2-two)*(five*c2 - one)*c2
     &            -((125.d0/6.d0))*c1*(five*c2 - two)*(five*c2 - one)*c2
              shl(2,18,lb)= ((625.d0/6.d0)*(five*c2 - one))*c1*c2*c3
     &              +((625.d0/6.d0)*(five*c2 - two))*c1*c2*c3
     &           +((125.d0/6.d0)*(five*c2 - two))*(five*c2 - one)*c1*c3
     &           -((125.d0/6.d0))*c1*(five*c2 - two)*(five*c2 - one)*c2
        shl(3,18,lb)=((125.d0/6.d0)*(five*c2-two))*(five*c2-one)*c1*c2*c3
              
              shl(1,19,lb)= -((625.d0/4.d0)*(five*c2 - one))*c1*c2*c3
     &           +((125.d0/4.d0)*(five*c3 - one))*c3*(five*c2 - one)*c2
     &           -((125.d0/4.d0)*(five*c2 - one))*(five*c3 - one)*c1*c2
              shl(2,19,lb)= ((625.d0/4.d0)*(five*c3 - one))*c1*c2*c3
     &              -((625.d0/4.d0)*(five*c2 - one))*c1*c2*c3
     &           +((125.d0/4.d0)*(five*c2 - one))*(five*c3 - one)*c1*c3
     &           -((125.d0/4.d0)*(five*c2 - one))*(five*c3 - one)*c1*c2
      shl(3,19,lb)= ((125.d0/4.d0)*(five*c2-one))*(five*c3-one)*c1*c2*c3
              
              shl(1,20,lb)= -((625.d0/6.d0)*(five*c3 - one))*c1*c2*c3
     &              -((625.d0/6.d0)*(five*c3 - two))*c1*c2*c3
     &           +((125.d0/6.d0)*(five*c3 - two))*(five*c3 - one)*c3*c2
     &           -((125.d0/6.d0)*(five*c3 - two))*(five*c3 - one)*c1*c2
              shl(2,20,lb)= -((625.d0/6.d0)*(five*c3 - one))*c1*c2*c3
     &              -((625.d0/6.d0)*(five*c3 - two))*c1*c2*c3
     &           +((125.d0/6.d0)*(five*c3 - two))*(five*c3 - one)*c3*c1
     &           -((125.d0/6.d0)*(five*c3 - two))*(five*c3 - one)*c1*c2
      shl(3,20,lb)= ((125.d0/6.d0)*(five*c3-two))*(five*c3-one)*c1*c2*c3

              shl(1,21,lb)= ((625.d0/4.d0)*(five*c3 - one))*c1*c2*c3
     &              -((625.d0/4.d0)*(five*c1 - one))*c1*c2*c3
     &           +((125.d0/4.d0)*(five*c1 - one))*(five*c3 - one)*c2*c3
     &           -((125.d0/4.d0)*(five*c1 - one))*(five*c3 - one)*c1*c2
              shl(2,21,lb)= -((625.d0/4.d0)*(five*c1 - one))*c1*c2*c3
     &           +((125.d0/4.d0)*(five*c3 - one))*c3*(five*c1 - one)*c1
     &           -((125.d0/4.d0)*(five*c1 - one))*(five*c3 - one)*c1*c2
      shl(3,21,lb)= ((125.d0/4.d0)*(five*c1-one))*(five*c3-one)*c1*c2*c3

            end if
          
  200     continue 
  100     continue
c
c 
           return 
           end
c**** new ********************************************************************** 
      subroutine shlq(shl,w,nint,nen)
c
c.... program to calculate integration-rule weights, shape functions
c        and local derivatives for a four-node quadrilateral element
c
c               s,t = local element coordinates ("xi", "eta", resp.)       
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c              w(l) = integration-rule weight
c                 i = local node number
c                 l = integration point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension wone(8),raone(8)
      dimension shl(3,nen,*),w(*),ra(64),sa(64)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      data  five9/0.5555555555555555d0/,eight9/0.8888888888888888d0/
      data r1/0.d00/,w1/2.d00/,
     &     r2/0.577350269189626d00/,w2/1.d00/,
     &     r3a/0.774596669241483d00/,w3a/0.555555555555556d00/,
     &     r3b/0.d00/,w3b/0.888888888888889d00/,
     &     r4a/0.861136311594053d00/,w4a/0.347854845137454d00/,
     &     r4b/0.339981043584856d00/,w4b/0.652145154862546d00/
c
      if (nint.eq.1) then
         wone(1)  = two
         raone(1) = zero
	nintx=1
	ninty=1
      endif
c
c
      if (nint.eq.4) then
         wone(1) = one
         wone(2) = one
         raone(1)=-.577350269189626
         raone(2)= .577350269189625
	nintx=2
	ninty=2
      endif
c
c
      if (nint.eq.9) then
         wone(1) = five9
         wone(2) = five9
         wone(3) = eight9
         raone(1)=-.774596669241483
         raone(2)= .774596669241483
         raone(3)= zero
	nintx=3
	ninty=3
      endif
c
      if (nint.eq.16) then
         wone(1) = .347854845137454
         wone(2) = .347854845137454
         wone(3) = .652145154862546
         wone(4) = .652145154862546
         raone(1)=-.861136311594053
         raone(2)= .861136311594053
         raone(3)=-.339981043584856
         raone(4)= .339981043584856
	nintx=4
	ninty=4
      endif
c
c
       if(nint.eq.25) then
        wone(1) = .236926885056189
        wone(2) = .236926885056189
        wone(3) = .478628670499366
        wone(4) = .478628670499366
        wone(5) = .568888888888888
        raone(1)=-.906179845938664
        raone(2)= .906179845938664
        raone(3)=-.538469310105683
        raone(4)= .538469310105683
        raone(5)= zero
	nintx=5
	ninty=5
       endif
c
       if(nint.eq.36) then
         wone(1) = .171324492397170
         wone(2) = .171324492397170
         wone(3) = .360761573048139
         wone(4) = .360761573048139
         wone(5) = .467913934572691
         wone(6) = .467913934572691
         raone(1)=-.932469514203152
         raone(2)= .932469514203152
         raone(3)=-.661209386466265
         raone(4)= .661209386466365
         raone(5)=-.238619186083197
         raone(6)= .238619186083197
	nintx=6
	ninty=6
       endif
c
c
       if(nint.eq.49) then
         wone(1) = .129484966168870
         wone(2) = .129484966168870
         wone(3) = .279705391489277 
         wone(4) = .279705391489277
         wone(5) = .381830050505119
         wone(6) = .381830050505119
         wone(7) = .417959183673469
         raone(1)=-.949107912342759
         raone(2)= .949107912342759
         raone(3)=-.741531185599394
         raone(4)= .741531185599394
         raone(5)=-.405845151377397
         raone(6)= .405845151377397
         raone(7)= zero
	nintx=7
	ninty=7
       endif
c
c
       if(nint.eq.64) then
         wone(1) = .101228536290376
         wone(2) = .101228536290376
         wone(3) = .222381034453374
         wone(4) = .222381034453374
         wone(5) = .313706645877887
         wone(6) = .313706645877887
         wone(7) = .362683783378362
         wone(8) = .362683783378362
         raone(1)=-.960289856497536
         raone(2)= .960289856497536
         raone(3)=-.796666477413627
         raone(4)= .796666477413627
         raone(5)=-.525532409916329
         raone(6)= .525532409916329
         raone(7)=-.183434642495650
         raone(8)= .183434642495650
	nintx=8
	ninty=8
       endif
c
c
       l=0
       do ly=1,ninty
	  do lx=1,nintx
	    l = l+1
	      w(l) = wone(lx)*wone(ly)
	      ra(l) = raone(lx)
	      sa(l) = raone(ly)
          end do
       end do
c
       do 200 l=1,nint
c
            r=ra(l)
            s=sa(l)
c
c
        if(nen.eq.1) then
            shl(1,1,l)=zero
            shl(2,1,l)=zero
            shl(3,1,l)=one
        end if
c
	if(nen.eq.4) then
            f1 = pt5*(one-r)
            f2 = pt5*(one+r)
	      fx1=-pt5
	      fx2= pt5
            g1 = pt5*(one-s)
            g2 = pt5*(one+s)
	      gx1=-pt5
	      gx2= pt5
            shl(1,1,l)=fx1*g1
            shl(2,1,l)=f1*gx1
            shl(3,1,l)=f1*g1
            shl(1,2,l)=fx2*g1
            shl(2,2,l)=f2*gx1
            shl(3,2,l)=f2*g1
            shl(1,3,l)=fx2*g2
            shl(2,3,l)=f2*gx2
            shl(3,3,l)=f2*g2
            shl(1,4,l)=fx1*g2
            shl(2,4,l)=f1*gx2
            shl(3,4,l)=f1*g2
	end if
c
c
        if(nen.eq.9) then
            f1 = -pt5*(one-r)*r
            f2 =  pt5*(one+r)*r
	      f3 = (one+r)*(one-r)
c
c	      
	      f1x = pt5*r - pt5*(one-r)
	      f2x = pt5*r + pt5*(one+r)
	      f3x =-two*r

            g1 =-pt5*(one-s)*s
            g2 = pt5*(one+s)*s
	      g3 = (one+s)*(one-s)
c
	      g1x = pt5*s - pt5*(one-s)
	      g2x = pt5*s + pt5*(one+s)
	      g3x =-two*s
c
c
            shl(3,1,l)=f1*g1
            shl(3,2,l)=f2*g1
            shl(3,3,l)=f2*g2
            shl(3,4,l)=f1*g2
c
                  shl(3,5,l)=f3*g1
                  shl(3,6,l)=f2*g3
                  shl(3,7,l)=f3*g2
                  shl(3,8,l)=f1*g3
                  shl(3,9,l)=f3*g3
c
            shl(1,1,l)=f1x*g1
            shl(1,2,l)=f2x*g1
            shl(1,3,l)=f2x*g2
            shl(1,4,l)=f1x*g2
c
                  shl(1,5,l)=f3x*g1
                  shl(1,6,l)=f2x*g3
                  shl(1,7,l)=f3x*g2
                  shl(1,8,l)=f1x*g3
                  shl(1,9,l)=f3x*g3
c
            shl(2,1,l)=f1*g1x
            shl(2,2,l)=f2*g1x
            shl(2,3,l)=f2*g2x
            shl(2,4,l)=f1*g2x
c
                  shl(2,5,l)=f3*g1x
                  shl(2,6,l)=f2*g3x
                  shl(2,7,l)=f3*g2x
                  shl(2,8,l)=f1*g3x
                  shl(2,9,l)=f3*g3x
c
        end if
c
c
        if (nen.eq.16) then
                  onemrsq=one-r*r
                  onemssq=one-s*s
                  onep3r=one+three*r
                  onem3r=one-three*r
                  onep3s=one+three*s
                  onem3s=one-three*s
		  f1=-1.d00/16.d00*(9.d00*r*r-1.d00)*(r-1.d00)
		  f2=9.d00/16.d00*(1.d00-r*r)*onem3r
		  f3=9.d00/16.d00*(1.d00-r*r)*onep3r
		  f4=1.d00/16.d00*(9.d00*r*r-1.d00)*(r+1.d00)
c
		  f1x=-1.d00/16.d00*(18.d00*r)*(r-1.d00)
     &               -1.d00/16.d00*(9.d00*r*r-1.d00)
		  f2x=9.d00/16.d00*(-2.d00*r)*onem3r
     &               -9.d00/16.d00*(1.d00-r*r)*3.d00
		  f3x=9.d00/16.d00*(-2.d00*r)*onep3r
     &               -9.d00/16.d00*(r*r-1.d00)*3.d00
      	  f4x=1.d00/16.d00*(18.d00*r)*(r+1.d00)
     &		       +1.d00/16.d00*(9.d00*r*r-1.d00)
c
            g1=-1.d00/16.d00*(9.d00*s*s-1.d00)*(s-1.d00)
		  g2=9.d00/16.d00*(1.d00-s*s)*onem3s
		  g3=9.d00/16.d00*(1.d00-s*s)*onep3s
		  g4=1.d00/16.d00*(9.d00*s*s-1.d00)*(s+1.d00)
c
		  g1x=-1.d00/16.d00*(18.d00*s)*(s-1.d00)
     &               -1.d00/16.d00*(9.d00*s*s-1.d00)
		  g2x=9.d00/16.d00*(-2.d00*s)*onem3s
     &               -9.d00/16.d00*(1.d00-s*s)*3.d00
		  g3x=9.d00/16.d00*(-2.d00*s)*onep3s
     &               -9.d00/16.d00*(s*s-1.d00)*3.d00
            g4x=1.d00/16.d00*(18.d00*s)*(s+1.d00)
     &		       +1.d00/16.d00*(9.d00*s*s-1.d00)
c
           shl(3,1,l)=f1*g1
	     shl(3,2,l)=f4*g1
	     shl(3,3,l)=f4*g4
	     shl(3,4,l)=f1*g4
c	     
           shl(3,5,l)=f2*g1
	     shl(3,6,l)=f3*g1
	     shl(3,7,l)=f4*g2
	     shl(3,8,l)=f4*g3
	     shl(3,9,l)=f3*g4
	     shl(3,10,l)=f2*g4
	     shl(3,11,l)=f1*g3
	     shl(3,12,l)=f1*g2
	     shl(3,13,l)=f2*g2
	     shl(3,14,l)=f3*g2
	     shl(3,15,l)=f3*g3
	     shl(3,16,l)=f2*g3
c
c
           shl(1,1,l)=f1x*g1
	     shl(1,2,l)=f4x*g1
	     shl(1,3,l)=f4x*g4
	     shl(1,4,l)=f1x*g4
c	     
           shl(1,5,l)=f2x*g1
	     shl(1,6,l)=f3x*g1
	     shl(1,7,l)=f4x*g2
	     shl(1,8,l)=f4x*g3
	     shl(1,9,l)=f3x*g4
	     shl(1,10,l)=f2x*g4
	     shl(1,11,l)=f1x*g3
	     shl(1,12,l)=f1x*g2
	     shl(1,13,l)=f2x*g2
	     shl(1,14,l)=f3x*g2
	     shl(1,15,l)=f3x*g3
	     shl(1,16,l)=f2x*g3
c
c
           shl(2,1,l)=f1*g1x
	     shl(2,2,l)=f4*g1x
	     shl(2,3,l)=f4*g4x
	     shl(2,4,l)=f1*g4x
c	     
           shl(2,5,l)=f2*g1x
	     shl(2,6,l)=f3*g1x
	     shl(2,7,l)=f4*g2x
	     shl(2,8,l)=f4*g3x
	     shl(2,9,l)=f3*g4x
	     shl(2,10,l)=f2*g4x
	     shl(2,11,l)=f1*g3x
	     shl(2,12,l)=f1*g2x
	     shl(2,13,l)=f2*g2x
	     shl(2,14,l)=f3*g2x
	     shl(2,15,l)=f3*g3x
	     shl(2,16,l)=f2*g3x
c
        end if
  200  continue
c
      return
      end
c**** new c**** new ********************************************************************** 
      subroutine shlqpk(shl,w,nint,nen)
c
c.... program to calculate integration-rule weights, shape functions
c        and local derivatives for a four-node quadrilateral element
c
c               s,t = local element coordinates ("xi", "eta", resp.)       
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c              w(l) = integration-rule weight
c                 i = local node number
c                 l = integration point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension wone(8),raone(8),xaone(8)
	dimension shlone(2,8,8),inod(8,8)
      dimension shl(3,nen,*),w(*),ra(64),sa(64)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      data  five9/0.5555555555555555d0/,eight9/0.8888888888888888d0/
      data r1/0.d00/,w1/2.d00/,
     &     r2/0.577350269189626d00/,w2/1.d00/,
     &     r3a/0.774596669241483d00/,w3a/0.555555555555556d00/,
     &     r3b/0.d00/,w3b/0.888888888888889d00/,
     &     r4a/0.861136311594053d00/,w4a/0.347854845137454d00/,
     &     r4b/0.339981043584856d00/,w4b/0.652145154862546d00/
c
      if (nint.eq.1) then
         wone(1)  = two
         raone(1) = zero
	nintx=1
	ninty=1
      endif
c
      if (nen.eq.1) then
         xaone(1) = zero
	 nenx=1
	 neny=1
      end if
c
      if (nint.eq.4) then
         wone(1) = one
         wone(2) = one
         raone(1)=-.577350269189626
         raone(2)= .577350269189625
	nintx=2
	ninty=2
      endif
c
      if (nen.eq.4) then
         xaone(1) = -one
         xaone(2) =  one
	nenx=2
	neny=2
      endif
c
c
c
      if (nint.eq.9) then
         wone(1) = five9
         wone(2) = five9
         wone(3) = eight9
         raone(1)=-.774596669241483
         raone(2)= .774596669241483
         raone(3)= zero
	nintx=3
	ninty=3
      endif
c 
      if(nen.eq.9) then
         xaone(1)= -one
         xaone(2)= one
         xaone(3)= zero
	nenx=3
	neny=3
      endif
c
      if (nint.eq.16) then
         wone(1) = .347854845137454
         wone(2) = .347854845137454
         wone(3) = .652145154862546
         wone(4) = .652145154862546
         raone(1)=-.861136311594053
         raone(2)= .861136311594053
         raone(3)=-.339981043584856
         raone(4)= .339981043584856
	nintx=4
	ninty=4
      endif
c
      if (nen.eq.16) then
         xaone(1) = -one
         xaone(2) = one
         xaone(3) = -.333333333333333
         xaone(4) =  .333333333333333
	nenx=4
	neny=4
         endif
c
c
c
       if(nint.eq.25) then
        wone(1) = .236926885056189
        wone(2) = .236926885056189
        wone(3) = .478628670499366
        wone(4) = .478628670499366
        wone(5) = .568888888888888
        raone(1)=-.906179845938664
        raone(2)= .906179845938664
        raone(3)=-.538469310105683
        raone(4)= .538469310105683
        raone(5)= zero
	nintx=5
	ninty=5
       endif
c
       if(nen.eq.25) then
         xaone(1)= -one 
         xaone(2)=  one 
         xaone(3)= -pt5
         xaone(4)= zero
         xaone(5)= pt5
	nenx=5
	neny=5
       endif
c
       if(nint.eq.36) then
         wone(1) = .171324492397170
         wone(2) = .171324492397170
         wone(3) = .360761573048139
         wone(4) = .360761573048139
         wone(5) = .467913934572691
         wone(6) = .467913934572691
         raone(1)=-.932469514203152
         raone(2)= .932469514203152
         raone(3)=-.661209386466265
         raone(4)= .661209386466365
         raone(5)=-.238619186083197
         raone(6)= .238619186083197
	nintx=6
	ninty=6
        endif
c
        if(nen.eq.36) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.600000000000000
         xaone(4) = -.200000000000000
         xaone(5) =  .200000000000000
         xaone(6) =  .600000000000000
	nenx=6
	neny=6
        endif
c
       if(nint.eq.49) then
         wone(1) = .129484966168870
         wone(2) = .129484966168870
         wone(3) = .279705391489277 
         wone(4) = .279705391489277
         wone(5) = .381830050505119
         wone(6) = .381830050505119
         wone(7) = .417959183673469
         raone(1)=-.949107912342759
         raone(2)= .949107912342759
         raone(3)=-.741531185599394
         raone(4)= .741531185599394
         raone(5)=-.405845151377397
         raone(6)= .405845151377397
         raone(7)= zero
	nintx=7
	ninty=7
        endif
c
        if(nen.eq.49) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.666666666666666
         xaone(4) = -.333333333333333
         xaone(5) = zero
         xaone(6) =  .333333333333333
         xaone(7) =  .666666666666666
	nenx=7
	neny=7
        endif
c
       if(nint.eq.64) then
         wone(1) = .101228536290376
         wone(2) = .101228536290376
         wone(3) = .222381034453374
         wone(4) = .222381034453374
         wone(5) = .313706645877887
         wone(6) = .313706645877887
         wone(7) = .362683783378362
         wone(8) = .362683783378362
         raone(1)=-.960289856497536
         raone(2)= .960289856497536
         raone(3)=-.796666477413627
         raone(4)= .796666477413627
         raone(5)=-.525532409916329
         raone(6)= .525532409916329
         raone(7)=-.183434642495650
         raone(8)= .183434642495650
	nintx=8
	ninty=8
        endif
        if(nen.eq.64) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -0.71428571428571
         xaone(4) = -0.42857142857143
         xaone(5) = -0.14285714285714
         xaone(6) = 0.14285714285714
         xaone(7) = 0.42857142857143
         xaone(8) = 0.71428571428571
	nenx=8
	neny=8
        endif
c
c
	if(nen.eq.1) inod(1,1) = 1
	if(nen.eq.4) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(1,2) = 4
	  inod(2,2) = 3
	end if
c
	if(nen.eq.9) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 7
c
	  inod(1,3) = 8
	  inod(2,3) = 6
	  inod(3,3) = 9
	end if
c
c
c
	if(nen.eq.16) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 10
	  inod(4,2) = 9
c
	  inod(1,3) = 12
	  inod(2,3) = 7
	  inod(3,3) = 13
	  inod(4,3) = 14
c
	  inod(1,4) = 11
	  inod(2,4) = 8
	  inod(3,4) = 16
	  inod(4,4) = 15
	end if
c
c
	if(nen.eq.25) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 13
	  inod(4,2) = 12
	  inod(5,2) = 11
c
	  inod(1,3) = 16
	  inod(2,3) = 8
	  inod(3,3) = 17
	  inod(4,3) = 18
	  inod(5,3) = 19
c
	  inod(1,4) = 15
	  inod(2,4) = 9
	  inod(3,4) = 24
	  inod(4,4) = 25
	  inod(5,4) = 20
c	  
	  inod(1,5) = 14
	  inod(2,5) = 10
	  inod(3,5) = 23
	  inod(4,5) = 22
	  inod(5,5) = 21
	end if
c
c
	if(nen.eq.36) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 16
	  inod(4,2) = 15
	  inod(5,2) = 14
	  inod(6,2) = 13
c
	  inod(1,3) = 20
	  inod(2,3) = 9
	  inod(3,3) = 21
	  inod(4,3) = 22
	  inod(5,3) = 23
	  inod(6,3) = 24
c
	  inod(1,4) = 19
	  inod(2,4) = 10
	  inod(3,4) = 32
	  inod(4,4) = 33
	  inod(5,4) = 34
	  inod(6,4) = 25
c	  
	  inod(1,5) = 18
	  inod(2,5) = 11
	  inod(3,5) = 31
	  inod(4,5) = 36
	  inod(5,5) = 35
	  inod(6,5) = 26
c
	  inod(1,6) = 17
	  inod(2,6) = 12
	  inod(3,6) = 30
	  inod(4,6) = 29
	  inod(5,6) = 28
	  inod(6,6) = 27
c	  
	end if
c
c
	if(nen.eq.49) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
	  inod(7,1) = 9
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 19
	  inod(4,2) = 18
	  inod(5,2) = 17
	  inod(6,2) = 16
	  inod(7,2) = 15
c
	  inod(1,3) = 24
	  inod(2,3) = 10
	  inod(3,3) = 25
	  inod(4,3) = 26
	  inod(5,3) = 27
	  inod(6,3) = 28
	  inod(7,3) = 29
c
	  inod(1,4) = 23
	  inod(2,4) = 11
	  inod(3,4) = 40
	  inod(4,4) = 41
	  inod(5,4) = 42
	  inod(6,4) = 43
	  inod(7,4) = 30
c	  
	  inod(1,5) = 22
	  inod(2,5) = 12
	  inod(3,5) = 39
	  inod(4,5) = 48
	  inod(5,5) = 49
	  inod(6,5) = 44
	  inod(7,5) = 31
c
	  inod(1,6) = 21
	  inod(2,6) = 13
	  inod(3,6) = 38
	  inod(4,6) = 47
	  inod(5,6) = 46
	  inod(6,6) = 45
	  inod(7,6) = 32
c
	  inod(1,7) = 20
	  inod(2,7) = 14
	  inod(3,7) = 37
	  inod(4,7) = 36
	  inod(5,7) = 35
	  inod(6,7) = 34
	  inod(7,7) = 33
c	  
	end if
c
c
	if(nen.eq.64) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
	  inod(7,1) = 9
	  inod(8,1) = 10
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 22
	  inod(4,2) = 21
	  inod(5,2) = 20
	  inod(6,2) = 19
	  inod(7,2) = 18
	  inod(8,2) = 17
c
	  inod(1,3) = 28
	  inod(2,3) = 11
	  inod(3,3) = 29
	  inod(4,3) = 30
	  inod(5,3) = 31
	  inod(6,3) = 32
	  inod(7,3) = 33
	  inod(8,3) = 34
c
	  inod(1,4) = 27
	  inod(2,4) = 12
	  inod(3,4) = 48
	  inod(4,4) = 49
	  inod(5,4) = 50
	  inod(6,4) = 51
	  inod(7,4) = 52
	  inod(8,4) = 35
c	  
	  inod(1,5) = 26
	  inod(2,5) = 13
	  inod(3,5) = 47
	  inod(4,5) = 60
	  inod(5,5) = 61
	  inod(6,5) = 62
	  inod(7,5) = 53
	  inod(8,5) = 36
c
	  inod(1,6) = 25
	  inod(2,6) = 14
	  inod(3,6) = 46
	  inod(4,6) = 59
	  inod(5,6) = 64
	  inod(6,6) = 63
	  inod(7,6) = 54
	  inod(8,6) = 37
c
	  inod(1,7) = 24
	  inod(2,7) = 15
	  inod(3,7) = 45
	  inod(4,7) = 58
	  inod(5,7) = 57
	  inod(6,7) = 56
	  inod(7,7) = 55
	  inod(8,7) = 38
c
	  inod(1,8) = 23
	  inod(2,8) = 16
	  inod(3,8) = 44
	  inod(4,8) = 43
	  inod(5,8) = 42
	  inod(6,8) = 41
	  inod(7,8) = 40
	  inod(8,8) = 39
c	  
	end if
c
c
      do 100 l = 1, nintx
         r = raone(l)
c
        if(nenx.eq.1) then
        shlone(1,1,l) = zero
        shlone(2,1,l) = one
        go to 100
        endif
c
        do 50 i = 1, nenx
         aa = one
         bb = one
         aax = zero
         do 40 j =1, nenx
          daj = one
          if (i .ne. j)then
          aa = aa * ( r - xaone(j))
          bb = bb * ( xaone(i) - xaone(j))
          do 30 k = 1, nenx
           if(k .ne. i .and. k .ne. j) daj = daj * ( r - xaone(k))
   30     continue
          aax =aax + daj
          endif
   40    continue
        shlone(2,i,l) = aa/bb
        shlone(1,i,l) = aax/bb
   50  continue
c
  100  continue      
c
         l=0
	   do ly=1,ninty
	   do lx=1,nintx
	   l = l+1
	      w(l) = wone(lx)*wone(ly)
         end do
	   end do
c
      l=0
	do ly=1,ninty
      do lx=1,nintx
	  l = l+1
	  r=ra(l)
	  s=sa(l)
	 do iy=1,neny
	 do ix=1,nenx
         j = inod(ix,iy)
	   shl(1,j,l) = shlone(1,ix,lx)*shlone(2,iy,ly)
	   shl(2,j,l) = shlone(2,ix,lx)*shlone(1,iy,ly)
	   shl(3,j,l) = shlone(2,ix,lx)*shlone(2,iy,ly)
	 end do
	 end do
	end do
	end do
c
      return
      end
c**** new ********************************************************************** 
      subroutine shlqpbk(shl,nen,nside,nnods,nints)
c
c.... program to calculate integration-rule weights, shape functions
c        and local derivatives for a four-node quadrilateral element
c
c               s,t = local element coordinates ("xi", "eta", resp.)       
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c              w(l) = integration-rule weight
c                 i = local node number
c                 l = integration point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension raone(8),xaone(8),paone(8)
	dimension shlx(2,8),shly(2,8),inod(8,8)
      dimension shl(3,nen,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      data  five9/0.5555555555555555d0/,eight9/0.8888888888888888d0/
      data r1/0.d00/,w1/2.d00/,
     &     r2/0.577350269189626d00/,w2/1.d00/,
     &     r3a/0.774596669241483d00/,w3a/0.555555555555556d00/,
     &     r3b/0.d00/,w3b/0.888888888888889d00/,
     &     r4a/0.861136311594053d00/,w4a/0.347854845137454d00/,
     &     r4b/0.339981043584856d00/,w4b/0.652145154862546d00/
c

      if (nints.eq.1) then
         raone(1) = zero
         
         paone(1) = zero
      endif
c
      if(nnods.eq.1) then
        xaone(1) = zero
      end if
c
      if (nints.eq.2) then
         raone(1)=-.577350269189626
         raone(2)= .577350269189625
c
         paone(1)= .577350269189626
         paone(2)=-.577350269189625
c
      endif
c
      if (nnods.eq.2) then
         xaone(1) = -one
         xaone(2) =  one
      endif
c
c
      if (nints.eq.3) then
         raone(1)=-.774596669241483
         raone(2)= zero
         raone(3)= .774596669241483
c
         paone(1)= .774596669241483
         paone(2)= zero
         paone(3)=-.774596669241483
c
      endif
c 
      if(nnods.eq.3) then
         xaone(1)= -one
         xaone(2)= one
         xaone(3)= zero
      endif
c
      if (nints.eq.4) then
         raone(1)=-.861136311594053
         raone(2)=-.339981043584856
         raone(3)= .339981043584856
         raone(4)= .861136311594053
c
         paone(1)= .861136311594053
         paone(2)= .339981043584856
         paone(3)=-.339981043584856
         paone(4)=-.861136311594053
c
      endif
c
      if (nnods.eq.4) then
         xaone(1) = -one
         xaone(2) = one
         xaone(3) = -.333333333333333
         xaone(4) =  .333333333333333
      endif
c
c
       if(nints.eq.5) then
        raone(1)=-.906179845938664
        raone(2)=-.538469310105683
        raone(3)= zero
        raone(4)= .538469310105683
        raone(5)= .906179845938664
c
        paone(1)= .906179845938664
        paone(2)= .538469310105683
        paone(3)= zero
        paone(4)=-.538469310105683
        paone(5)=-.906179845938664
c
       endif
c
       if(nnods.eq.5) then
         xaone(1)= -one 
         xaone(2)=  one 
         xaone(3)= -pt5
         xaone(4)= zero
         xaone(5)= pt5
       endif
c
       if(nints.eq.6) then
         raone(1)=-.932469514203152
         raone(2)=-.661209386466265
         raone(3)=-.238619186083197
         raone(4)= .238619186083197
         raone(5)= .661209386466365
         raone(6)= .932469514203152
c
         paone(1)= .932469514203152
         paone(2)= .661209386466265
         paone(3)= .238619186083197
         paone(4)=-.238619186083197
         paone(5)=-.661209386466365
         paone(6)=-.932469514203152
c
        endif
c
        if(nnods.eq.6) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.600000000000000
         xaone(4) = -.200000000000000
         xaone(5) =  .200000000000000
         xaone(6) =  .600000000000000
        endif
c
       if(nints.eq.7) then
         raone(1)=-.949107912342759
         raone(2)=-.741531185599394
         raone(3)=-.405845151377397
         raone(4)= zero
         raone(5)= .405845151377397
         raone(6)= .741531185599394
         raone(7)= .949107912342759
c
         paone(1)= .949107912342759
         paone(2)= .741531185599394
         paone(3)= .405845151377397
         paone(4)= zero
         paone(5)=-.405845151377397
         paone(6)=-0.741531185599394
         paone(7)=-.949107912342759
c
        endif
c
        if(nnods.eq.7) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.666666666666666
         xaone(4) = -.333333333333333
         xaone(5) = zero
         xaone(6) =  .333333333333333
         xaone(7) =  .666666666666666
        endif
c
       if(nints.eq.8) then
         raone(1)=-.960289856497536
         raone(2)=-.796666477413627
         raone(3)=-.525532409916329
         raone(4)=-.183434642495650
         raone(5)= .183434642495650
         raone(6)= .525532409916329
         raone(7)= .796666477413627
         raone(8)= .960289856497536
c
         paone(1)= .960289856497536
         paone(2)= .796666477413627
         paone(3)= .525532409916329
         paone(4)= .183434642495650
         paone(5)=-.183434642495650
         paone(6)=-.525532409916329
         paone(7)=-.796666477413627
         paone(8)=-.960289856497536
c
        endif
        if(nnods.eq.8) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -0.71428571428571
         xaone(4) = -0.42857142857143
         xaone(5) = -0.14285714285714
         xaone(6) = 0.14285714285714
         xaone(7) = 0.42857142857143
         xaone(8) = 0.71428571428571
        endif
c
c
	if(nen.eq.1) inod(1,1) = 1
	if(nen.eq.4) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(1,2) = 4
	  inod(2,2) = 3
	end if
c
	if(nen.eq.9) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 7
c
	  inod(1,3) = 8
	  inod(2,3) = 6
	  inod(3,3) = 9
	end if
c
	if(nen.eq.16) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 10
	  inod(4,2) = 9
c
	  inod(1,3) = 12
	  inod(2,3) = 7
	  inod(3,3) = 13
	  inod(4,3) = 14
c
	  inod(1,4) = 11
	  inod(2,4) = 8
	  inod(3,4) = 16
	  inod(4,4) = 15
	end if
c
c
	if(nen.eq.25) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 13
	  inod(4,2) = 12
	  inod(5,2) = 11
c
	  inod(1,3) = 16
	  inod(2,3) = 8
	  inod(3,3) = 17
	  inod(4,3) = 18
	  inod(5,3) = 19
c
	  inod(1,4) = 15
	  inod(2,4) = 9
	  inod(3,4) = 24
	  inod(4,4) = 25
	  inod(5,4) = 20
c	  
	  inod(1,5) = 14
	  inod(2,5) = 10
	  inod(3,5) = 23
	  inod(4,5) = 22
	  inod(5,5) = 21
	end if
c
	if(nen.eq.36) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 16
	  inod(4,2) = 15
	  inod(5,2) = 14
	  inod(6,2) = 13
c
	  inod(1,3) = 20
	  inod(2,3) = 9
	  inod(3,3) = 21
	  inod(4,3) = 22
	  inod(5,3) = 23
	  inod(6,3) = 24
c
	  inod(1,4) = 19
	  inod(2,4) = 10
	  inod(3,4) = 32
	  inod(4,4) = 33
	  inod(5,4) = 34
	  inod(6,4) = 25
c	  
	  inod(1,5) = 18
	  inod(2,5) = 11
	  inod(3,5) = 31
	  inod(4,5) = 36
	  inod(5,5) = 35
	  inod(6,5) = 26
c
	  inod(1,6) = 17
	  inod(2,6) = 12
	  inod(3,6) = 30
	  inod(4,6) = 29
	  inod(5,6) = 28
	  inod(6,6) = 27
c	  
	end if
c
	if(nen.eq.49) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
	  inod(7,1) = 9
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 19
	  inod(4,2) = 18
	  inod(5,2) = 17
	  inod(6,2) = 16
	  inod(7,2) = 15
c
	  inod(1,3) = 24
	  inod(2,3) = 10
	  inod(3,3) = 25
	  inod(4,3) = 26
	  inod(5,3) = 27
	  inod(6,3) = 28
	  inod(7,3) = 29
c
	  inod(1,4) = 23
	  inod(2,4) = 11
	  inod(3,4) = 40
	  inod(4,4) = 41
	  inod(5,4) = 42
	  inod(6,4) = 43
	  inod(7,4) = 30
c	  
	  inod(1,5) = 22
	  inod(2,5) = 12
	  inod(3,5) = 39
	  inod(4,5) = 48
	  inod(5,5) = 49
	  inod(6,5) = 44
	  inod(7,5) = 31
c
	  inod(1,6) = 21
	  inod(2,6) = 13
	  inod(3,6) = 38
	  inod(4,6) = 47
	  inod(5,6) = 46
	  inod(6,6) = 45
	  inod(7,6) = 32
c
	  inod(1,7) = 20
	  inod(2,7) = 14
	  inod(3,7) = 37
	  inod(4,7) = 36
	  inod(5,7) = 35
	  inod(6,7) = 34
	  inod(7,7) = 33
c	  
	end if
c
c
	if(nen.eq.64) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
	  inod(7,1) = 9
	  inod(8,1) = 10
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 22
	  inod(4,2) = 21
	  inod(5,2) = 20
	  inod(6,2) = 19
	  inod(7,2) = 18
	  inod(8,2) = 17
c
	  inod(1,3) = 28
	  inod(2,3) = 11
	  inod(3,3) = 29
	  inod(4,3) = 30
	  inod(5,3) = 31
	  inod(6,3) = 32
	  inod(7,3) = 33
	  inod(8,3) = 34
c
	  inod(1,4) = 27
	  inod(2,4) = 12
	  inod(3,4) = 48
	  inod(4,4) = 49
	  inod(5,4) = 50
	  inod(6,4) = 51
	  inod(7,4) = 52
	  inod(8,4) = 35
c	  
	  inod(1,5) = 26
	  inod(2,5) = 13
	  inod(3,5) = 47
	  inod(4,5) = 60
	  inod(5,5) = 61
	  inod(6,5) = 62
	  inod(7,5) = 53
	  inod(8,5) = 36
c
	  inod(1,6) = 25
	  inod(2,6) = 14
	  inod(3,6) = 46
	  inod(4,6) = 59
	  inod(5,6) = 64
	  inod(6,6) = 63
	  inod(7,6) = 54
	  inod(8,6) = 37
c
	  inod(1,7) = 24
	  inod(2,7) = 15
	  inod(3,7) = 45
	  inod(4,7) = 58
	  inod(5,7) = 57
	  inod(6,7) = 56
	  inod(7,7) = 55
	  inod(8,7) = 38
c
	  inod(1,8) = 23
	  inod(2,8) = 16
	  inod(3,8) = 44
	  inod(4,8) = 43
	  inod(5,8) = 42
	  inod(6,8) = 41
	  inod(7,8) = 40
	  inod(8,8) = 39
c	  
	end if
c
      lb=0
      do 200 ns=1,nside    
      do 200 l = 1, nints
      lb = lb + 1
      if(ns.eq.1) then
         r = raone(l)
         s = -one
      end if 
      if(ns.eq.2) then
         r = one
         s = raone(l)
      end if 
c
      if(ns.eq.3) then
         r = paone(l)
         s = one
      end if 
c
      if(ns.eq.4) then
         r = -one
         s = paone(l)
      end if 
c
        shlx(1,1) = zero
        shly(1,1) = zero
        shlx(2,1) = one
        shly(2,1) = one
      if(nnods.eq.1) go to 100
c
        do 50 i = 1, nnods
         aa = one
         bb = one
         cc = one
         aax= zero
         aay= zero
         do 40 j =1, nnods
          daj = one
          caj = one
          if (i .ne. j)then
          aa = aa * ( r - xaone(j))
          cc = cc * ( s - xaone(j))
          bb = bb * ( xaone(i) - xaone(j))
          do 30 k = 1, nnods
           if(k .ne. i .and. k .ne. j) daj = daj * ( r - xaone(k))
           if(k .ne. i .and. k .ne. j) caj = caj * ( s - xaone(k))
   30     continue
          aax = aax + daj
          aay = aay + caj
          endif
   40    continue
        shlx(1,i) = aax/bb
        shly(1,i) = aay/bb
        shlx(2,i) = aa/bb
        shly(2,i) = cc/bb
   50  continue
c
c
  100  continue 
c
	 do iy=1,nnods
	 do ix=1,nnods
         j = inod(ix,iy)
	   shl(1,j,lb) = shlx(1,ix)*shly(2,iy)
	   shl(2,j,lb) = shlx(2,ix)*shly(1,iy)
	   shl(3,j,lb) = shlx(2,ix)*shly(2,iy)
	 end do
	 end do
 200  continue
c
      return
      end
c*******************************************************************
      subroutine shgqs(xl,det,shl,shg,nint,nel,neg,quad,nen,
     &                 shlnode,nenode)
c
c.... program to calculate global derivatives of shape functions and
c        jacobian determinants for a  quadrilateral element
c
c        xl(j,i)    = global coordinates
c        det(l)     = jacobian determinant
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c        shg(1,i,l) = x-derivative of shape function
c        shg(2,i,l) = y-derivative of shape function
c        shg(3,i,l) = shl(3,i,l)
c        xs(i,j)    = jacobian matrix
c                 i = local node number or global coordinate number
c                 j = global coordinate number
c                 l = integration-point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      logical quad
      dimension xl(2,*),det(*),shl(3,nen,*),shg(3,nen,*),xs(2,2),
     &          shlnode(3,nenode,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c      
      call move(shg,shl,3*nen*nint)
c
      do 700 l=1,nint
c
      if (.not.quad) then
         do 100 i=1,3
         shg(i,3,l) = shl(i,3,l) + shl(i,4,l)
         shg(i,4,l) = zero
  100    continue
      endif
c
      do 300 j=1,2
      do 200 i=1,2
      xs(i,j) = rowdot(shlnode(i,1,l),xl(j,1),3,2,nenode)
  200 continue
  300 continue
c
      det(l) = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      if (det(l).le.zero) then
         write(iecho,1000) nel,neg
         stop
      endif
c
      do 500 j=1,2
      do 400 i=1,2
      xs(i,j) = xs(i,j)/det(l) 
  400 continue
  500 continue
c
      do 600 i=1,nen
        temp = xs(2,2)*shg(1,i,l) - xs(1,2)*shg(2,i,l)
        shg(2,i,l) = - xs(2,1)*shg(1,i,l) + xs(1,1)*shg(2,i,l)
        shg(1,i,l) = temp
  600 continue
c
  700 continue
c
      return
c
 1000 format(///,'shgqs - non-positive determinant - element ',i10,
     &          ' in element group  ',i10,'  shgqs')
      end
c**** new *********************************************************************
      subroutine shgtqsd(ien,xl,shl,shg,
     &                 nside,nint,nints,nel,neg,nen,
     &                 shlnode,nenode)
c
c.... program to calculate global derivatives of shape functions and
c        jacobian determinants for a  quadrilateral element
c
c        xl(j,i)    = global coordinates
c        det(l)     = jacobian determinant
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c        shg(1,i,l) = x-derivative of shape function
c        shg(2,i,l) = y-derivative of shape function
c        shg(3,i,l) = shl(3,i,l)
c        xs(i,j)    = jacobian matrix
c                 i = local node number or global coordinate number
c                 j = global coordinate number
c                 l = integration-point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension xl(2,*),shl(3,nen,*),shg(3,nen,*),xs(2,2),
     &          shlnode(3,nenode,*),ien(*)
      dimension igas(100)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi 
c
cc      call move(shg,shl,3*nen*nint)
c
c   renumeracao
c
      ns = 1
      nl1=ien(1)
      nl2=ien(2)
c      
      if(nl2.gt.nl1) then 
        do nn=1,nints
           ng = (ns-1)*nints + nn
           igas(ng) = ng
	  end do
	else
        do  nn=1,nints
            ng = (ns-1)*nints + nn
            ngs = (ns-1)*nints
	      igas(ng) = ngs + nints + 1 - nn
	    end do 
      end if 
c
c
      ns = 2
      nl1=ien(2)
      nl2=ien(3)
c      
      if(nl2.gt.nl1) then 
        do  nn=1,nints
           ng = (ns-1)*nints + nn
           igas(ng) = ng
	  end do
	else
        do  nn=1,nints
            ng = (ns-1)*nints + nn
            ngs = (ns-1)*nints
	      igas(ng) = ngs + nints + 1 - nn
	    end do 
      end if 
c
c
      if(nside.eq.3) then 
        ns = 3
        nl1=ien(3)
        nl2=ien(1)
c      
        if(nl2.gt.nl1) then 
          do  nn=1,nints
             ng = (ns-1)*nints + nn
             igas(ng) = ng
	    end do
	   else
           do  nn=1,nints
             ng = (ns-1)*nints + nn
             ngs = (ns-1)*nints
	       igas(ng) = ngs + nints + 1 - nn
	     end do 
         end if 
c
      end if 
c
c
      if(nside.eq.4) then 
        ns = 3
        nl1=ien(3)
        nl2=ien(4)
c      
        if(nl2.gt.nl1) then 
          do  nn=1,nints
             ng = (ns-1)*nints + nn
             igas(ng) = ng
	    end do
	  else
          do  nn=1,nints
              ng = (ns-1)*nints + nn
              ngs = (ns-1)*nints
	        igas(ng) = ngs + nints + 1 - nn
	    end do 
         end if 
c 
        ns = 4
        nl1=ien(4)
        nl2=ien(1)
c      
        if(nl2.gt.nl1) then 
          do  nn=1,nints
             ng = (ns-1)*nints + nn
             igas(ng) = ng
	    end do
	  else
          do  nn=1,nints
              ng = (ns-1)*nints + nn
              ngs = (ns-1)*nints
	        igas(ng) = ngs + nints + 1 - nn
	    end do 
         end if 
c 
      end if
c 
c
      do 700 l=1,nint
c
c
      do 300 j=1,2
      do 200 i=1,2
      xs(i,j) = rowdot(shlnode(i,1,l),xl(j,1),3,2,nenode)
  200 continue
  300 continue
c
      deta = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      if (deta.le.zero) then
         write(iecho,1000) nel,neg
         stop
      endif
c
      do 500 j=1,2
      do 400 i=1,2
      xs(i,j) = xs(i,j)/deta 
  400 continue
  500 continue
c
      do 600 i=1,nen
c
        lga = igas(l)
        shg(3,i,lga) = shl(3,i,l)        
        temp = xs(2,2)*shl(1,i,l) - xs(1,2)*shl(2,i,l)
        shg(2,i,lga) = - xs(2,1)*shl(1,i,l) + xs(1,1)*shl(2,i,l)
        shg(1,i,lga) = temp
  600 continue
c
  700 continue
c      
      return
c
 1000 format(///,'shgtqsd - non-positive determinant - element ',i10,
     &          ' in element group  ',i10,'  shgtqsd')
      end
c***********************************************************************      
      subroutine oneshgp(xl,det,shl,shlp,shgp,
     &                   nen,npars,nint,nesd,ns,nel,neg)
c
c.... program to calculate global derivatives of shape functions 
c        and jacobian determinants for the bi-dimensional,
c        elastic beam element
c
c           xl(j,l) = global coordinates of integration points
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = shape function
c        shg(1,i,l) = global ("arc-length") derivative of shape ftn
c        shg(2,i,l) = shl(2,i,l)
c            det(l) = euclidean length 
c                 i = local node number 
c                 j = global coordinate number
c                 l = integration-point number
c              nint = number of integration points
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension xl(nesd,*),det(*),shl(3,nen,*)
	dimension shlp(3,npars,*),shgp(3,npars,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
c
      do 400 l=1,nint
c
      det(l)=zero
      x1=0.d0
      x2=0.d0
      do 100 j=1,nen
      x1=x1+shl(1,j,l)*xl(1,j)
      x2=x2+shl(1,j,l)*xl(2,j)
100   continue
      det(l)=dsqrt(x1*x1+x2*x2)
c
      if (det(l).le.zero) then
         write(iecho,1000) ns,nel,neg
         stop
      endif
c
      do 300 i=1,npars
      shgp(1,i,l)=shlp(1,i,l)/det(l)
      shgp(2,i,l)=shlp(2,i,l)
300   continue
c
  400 continue
 1000 format(///,' oneshgp - non-positive determinant in side ',i10,/,
     &          ' in element ',i10,5x,'in oneshgp group  ',i10)
c
      return
      end
c**** new ********************************************************************** 
      subroutine shapesd(shl,nen,nenc)
c
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension raone(8),xaone(8)
      dimension shl(8,8)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      if (nen.eq.1) xaone(1) = zero
      if (nenc.eq.1) raone(1) = zero
c
c
      if (nen.eq.2) then
         xaone(1) = -one
         xaone(2) =  one
      endif
c
      if (nenc.eq.2) then
         raone(1) = -one
         raone(2) =  one
      endif
c
c 
      if(nen.eq.3) then
         xaone(1)= -one
         xaone(2)= one
         xaone(3)= zero
       endif
c 
      if(nenc.eq.3) then
         raone(1)= -one
         raone(2)= one
         raone(3)= zero
       endif
c
c
      if (nen.eq.4) then
         xaone(1) = -one
         xaone(2) = one
         xaone(3) = -.333333333333333
         xaone(4) =  .333333333333333
       endif
c
      if (nenc.eq.4) then
         raone(1) = -one
         raone(2) = one
         raone(3) = -.333333333333333
         raone(4) =  .333333333333333
       endif
c
c
       if(nen.eq.5) then
         xaone(1)= -one 
         xaone(2)=  one 
         xaone(3)= -pt5
         xaone(4)= zero
         xaone(5)= pt5
        endif
c
       if(nenc.eq.5) then
         raone(1)= -one 
         raone(2)=  one 
         raone(3)= -pt5
         raone(4)= zero
         raone(5)= pt5
       endif
c
c
        if(nen.eq.6) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.600000000000000
         xaone(4) = -.200000000000000
         xaone(5) =  .200000000000000
         xaone(6) =  .600000000000000
        endif
c
        if(nenc.eq.6) then
         raone(1) = -one
         raone(2) =  one
         raone(3) = -.600000000000000
         raone(4) = -.200000000000000
         raone(5) =  .200000000000000
         raone(6) =  .600000000000000
        endif
c
c
        if(nen.eq.7) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.666666666666666
         xaone(4) = -.333333333333333
         xaone(5) = zero
         xaone(6) =  .333333333333333
         xaone(7) =  .666666666666666
        endif
c
        if(nenc.eq.7) then
         raone(1) = -one
         raone(2) =  one
         raone(3) = -.666666666666666
         raone(4) = -.333333333333333
         raone(5) = zero
         raone(6) =  .333333333333333
         raone(7) =  .666666666666666
        endif
c
c
        if(nen.eq.8) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -0.71428571428571
         xaone(4) = -0.42857142857143
         xaone(5) = -0.14285714285714
         xaone(6) = 0.14285714285714
         xaone(7) = 0.42857142857143
         xaone(8) = 0.71428571428571
        endif
c
        if(nenc.eq.8) then
         raone(1) = -one
         raone(2) =  one
         raone(3) = -0.71428571428571
         raone(4) = -0.42857142857143
         raone(5) = -0.14285714285714
         raone(6) = 0.14285714285714
         raone(7) = 0.42857142857143
         raone(8) = 0.71428571428571
        endif
c
c
c
      do 100 l = 1, nenc
         r = raone(l)
c
        if(nen.eq.1) then
         shl(1,l) = zero
         go to 100
        endif
c
        do 50 i = 1, nen
         aa = one
         bb = one
         aax = zero
         do 40 j =1, nen
          daj = one
          if (i .ne. j)then
          aa = aa * ( r - xaone(j))
          bb = bb * ( xaone(i) - xaone(j))
          endif
   40    continue
        shl(i,l) = aa/bb
   50  continue
c
  100  continue      
c
      return
      end
c**** new ********************************************************************** 
      subroutine shapent(shl,nen,nenc)
c
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension shl(64,64),cl1(64),cl2(64),cl3(64)
      data   zero,pt1667,pt25,pt5 
     &      /0.0d0,0.1666666666666667d0,0.25d0,0.5d0/, 
     &       one,two,three,four,five,six 
     &      /1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0/ 
      data r1/0.33333333333333333333d00/ ,w1/1.d00/, 
     &     r2/0.5d00                   /,w2/0.3333333333333333333d00/
c

      pt1s3 = one/three
      pt2s3 = two/three
      pt1s4 = one/four
      pt2s4 = two/four
      pt3s4 = three/four
      pt1s5 = one/five
      pt2s5 = two/five
      pt3s5 = three/five
      pt4s5 = four/five
c      
      pt3s2 = three/two 
      pt9s2 = 9.d0/two
      pt27s2= 27.d0/two
      pt2s3 = two/three
      pt1s6 = one/six
      pt32s3= 32.d0/three
      pt8s3 = 8.d0/three
     
c
      if (nenc.eq.1) then 
            cl1(1)=r1 
            cl2(1)=r1 
            cl3(1)=one-r1-r1 
      end if 
c
      if(nenc.eq.3) then 
        cl1(1)=one 
        cl2(1)=zero 
        cl3(1)=zero
        
        cl1(2)=zero 
        cl2(2)=one 
        cl3(2)=zero
        
        cl1(3)=zero 
        cl2(3)=zero 
        cl3(3)=one 
      end if 
c 
      if(nenc.eq.6) then 
        cl1(1)=one 
        cl2(1)=zero 
        cl3(1)=zero
        
        cl1(2)=zero 
        cl2(2)=one 
        cl3(2)=zero
        
        cl1(3)=zero 
        cl2(3)=zero 
        cl3(3)=one 
        
        cl1(4)=r2 
        cl2(4)=r2 
        cl3(4)=zero 
        
        cl1(5)=zero 
        cl2(5)=r2 
        cl3(5)=r2
        
        cl1(6)=r2 
        cl2(6)=zero 
        cl3(6)=r2
        
      end if 
c
      if(nenc.eq.10) then
        cl1(1)=one 
        cl2(1)=zero 
        cl3(1)=zero
        
        cl1(2)=zero 
        cl2(2)=one 
        cl3(2)=zero
        
        cl1(3)=zero 
        cl2(3)=zero 
        cl3(3)=one 

        cl1(4)=pt2s3 
        cl2(4)=pt1s3 
        cl3(4)=zero
        
        cl1(5)=pt1s3 
        cl2(5)=pt2s3
        cl3(5)=zero 
        
        cl1(6)=zero 
        cl2(6)=pt2s3
        cl3(6)=pt1s3 
        
        cl1(7)=zero 
        cl2(7)=pt1s3
        cl3(7)=pt2s3 
        
        cl1(8)=pt1s3 
        cl2(8)=zero
        cl3(8)=pt2s3 
        
        cl1(9)=pt2s3 
        cl2(9)=zero
        cl3(9)=pt1s3 
        
        cl1(10)=r1 
        cl2(10)=r1
        cl3(10)=r1 

      end if
c
      if(nenc.eq.15) then
        cl1(1)=one 
        cl2(1)=zero 
        cl3(1)=zero
        
        cl1(2)=zero 
        cl2(2)=one
        cl3(2)=zero 
        
        cl1(3)=zero 
        cl2(3)=zero 
        cl3(3)=one 

        cl1(4)=pt3s4 
        cl2(4)=pt1s4 
        cl3(4)=zero
        
        cl1(5)=pt2s4 
        cl2(5)=pt2s4
        cl3(5)=zero 
        
        cl1(6)=pt1s4
        cl2(6)=pt3s4
        cl3(6)=zero 
        
        cl1(7)=zero 
        cl2(7)=pt3s4
        cl3(7)=pt1s4 
        
        cl1(8)=zero 
        cl2(8)=pt2s4
        cl3(8)=pt2s4 
        
        cl1(9)=zero 
        cl2(9)=pt1s4
        cl3(9)=pt3s4 
        
        cl1(10)=pt1s4 
        cl2(10)=zero
        cl3(10)=pt3s4

        cl1(11)=pt2s4 
        cl2(11)=zero
        cl3(11)=pt2s4

        cl1(12)=pt3s4 
        cl2(12)=zero
        cl3(12)=pt1s4

        cl1(13)=pt2s4 
        cl2(13)=pt1s4
        cl3(13)=pt1s4

        cl1(14)=pt1s4 
        cl2(14)=pt2s4
        cl3(14)=pt1s4

        cl1(15)=pt1s4 
        cl2(15)=pt1s4
        cl3(15)=pt2s4
      end if
c
      if(nenc.eq.21) then
        cl1(1)=one 
        cl2(1)=zero 
        cl3(1)=zero
        
        cl1(2)=zero 
        cl2(2)=one
        cl3(2)=zero 
        
        cl1(3)=zero 
        cl2(3)=zero 
        cl3(3)=one 

        cl1(4)=pt4s5 
        cl2(4)=pt1s5 
        cl3(4)=zero
        
        cl1(5)=pt3s5 
        cl2(5)=pt2s5
        cl3(5)=zero 
        
        cl1(6)=pt2s5
        cl2(6)=pt3s5
        cl3(6)=zero 
        
        cl1(7)=pt1s5
        cl2(7)=pt4s5
        cl3(7)=zero 
        
        cl1(8)=zero 
        cl2(8)=pt4s5
        cl3(8)=pt1s5 
        
        cl1(9)=zero 
        cl2(9)=pt3s5
        cl3(9)=pt2s5 
        
        cl1(10)=zero 
        cl2(10)=pt2s5
        cl3(10)=pt3s5

        cl1(11)=zero 
        cl2(11)=pt1s5
        cl3(11)=pt4s5

        cl1(12)=pt1s5 
        cl2(12)=zero
        cl3(12)=pt4s5

        cl1(13)=pt2s5 
        cl2(13)=zero
        cl3(13)=pt3s5

        cl1(14)=pt3s5 
        cl2(14)=zero
        cl3(14)=pt2s5

        cl1(15)=pt4s5 
        cl2(15)=zero
        cl3(15)=pt1s5
        
        cl1(16)=pt3s5 
        cl2(16)=pt1s5
        cl3(16)=pt1s5

        cl1(17)=pt2s5 
        cl2(17)=pt2s5
        cl3(17)=pt1s5

        cl1(18)=pt1s5 
        cl2(18)=pt3s5
        cl3(18)=pt1s5

        cl1(19)=pt1s5 
        cl2(19)=pt2s5
        cl3(19)=pt2s5

        cl1(20)=pt1s5 
        cl2(20)=pt1s5
        cl3(20)=pt3s5

        cl1(21)=pt2s5 
        cl2(21)=pt1s5
        cl3(21)=pt2s5    
      end if
c
      do 200 l=1,nenc 
c 
            c1 = cl1(l) 
            c2 = cl2(l) 
            c3 = cl3(l) 
!          acrescentei o if (2012-10-08)
!           interpolação linear (p=1 e nen=3)(2012-10-08)
            if(nen.eq.3) then
        shl(1,l)= c1
        shl(2,l)= c2
        shl(3,l)= c3
            end if 
!           interpolação quadrática (p=2 e nen=6)(2012-11-24)
            if(nen.eq.6) then
        shl(1,l)= (two*c1 - one)*c1
        shl(2,l)= (two*c2 - one)*c2
        shl(3,l)= (two*c3 - one)*c3
        shl(4,l)= four * c1 * c2
        shl(5,l)= four * c2 * c3
        shl(6,l)= four * c3 * c1
            end if
!           interpolação cúbica (p=3 e nen=10)(2012-11-24)
            if(nen.eq.10) then
           
        shl(1,l)= pt5*c1*(three*c1 - two)*(three*c1-one)
              
        shl(2,l)= pt5*c2*(three*c2 - two)*(three*c2-one)
              
        shl(3,l)= pt5*c3*(three*c3 - two)*(three*c3-one)  
              
        shl(4,l)= pt9s2*c1*c2*(three*c1-one)
              
        shl(5,l)= pt9s2*c1*c2*(three*c2-one)
              
        shl(6,l)= pt9s2*c2*c3*(three*c2-one)

        shl(7,l)= pt9s2*c2*c3*(three*c3-one)

        shl(8,l)= pt9s2*c3*c1*(three*c3-one)

        shl(9,l)= pt9s2*c3*c1*(three*c1-one)

        shl(10,l)= 27.d0*c1*c2*c3

            end if
!           interpolação quártica (p=4 e nen=15)(2012-11-24)
            if(nen.eq.15) then
           
        shl(1,l)= pt1s6*(four*c1-three)*(four*c1-two)*
     &                  (four*c1-one)*c1
              
        shl(2,l)= pt1s6*(four*c2-three)*(four*c2-two)*
     &                  (four*c2-one)*c2
              
        shl(3,l)= pt1s6*(four*c3-three)*(four*c3-two)*
     &                   (four*c3-one)*c3
              
        shl(4,l)= pt8s3*c2*(four*c1-two)*(four*c1-one)*c1
              
        shl(5,l)= four*c2*(four*c2-one)*(four*c1-one)*c1
              
        shl(6,l)= pt8s3*c1*(four*c2-two)*(four*c2-one)*c2

        shl(7,l)= pt8s3*c3*(four*c2-two)*(four*c2-one)*c2

        shl(8,l)= four*(four*c2-one)*c2*(four*c3-one)*c3

        shl(9,l)= pt8s3*c2*(four*c3-two)*(four*c3-one)*c3

        shl(10,l)= pt8s3*c1*(four*c3-two)*(four*c3-one)*c3

        shl(11,l)= 4.d0*(four*c3-one)*c3*(four*c1-one)*c1
              
        shl(12,l)= pt8s3*c3*(four*c1-two)*(four*c1-one)*c1
              
        shl(13,l)= 32.d0*c1*c2*c3*(four*c1-one)
              
        shl(14,l)= 32.d0*c1*c2*c3*(four*c2-one)

        shl(15,l)= 32.d0*c1*c2*c3*(four*c3-one)

            end if
!           interpolação quíntica (p=5 e nen=21)(2012-11-24)
            if(nen.eq.21) then
              
        shl(1,l)= (one/24.d0)*(five*c1-four)*(five*c1-three)
     &                    *(five*c1-two)*(five*c1-one)*c1
              
        shl(2,l)= (one/24.d0)*(five*c2-four)*(five*c2-three)
     &                    *(five*c2-two)*(five*c2-one)*c2
              
        shl(3,l)= (one/24.d0)*(five*c3-four)*(five*c3-three)
     &                    *(five*c3-two)*(five*c3-one)*c3
              
        shl(4,l)= (25.d0/24.d0)*(five*c1-three)
     &                    *(five*c1-two)*(five*c1-one)*c1*c2
              
        shl(5,l)= (25.d0/12.d0)*(five*c1-two)
     &                    *(five*c1-one)*c1*(five*c2-one)*c2
              
        shl(6,l)= (25.d0/12.d0)*(five*c1-one)
     &                    *c1*(five*c2-two)*(five*c2-one)*c2

        shl(7,l)= (25.d0/24.d0)*c1*
     &                (five*c2-three)*(five*c2-two)*(five*c2-one)*c2

        shl(8,l)= (25.d0/24.d0)*c3*
     &                 (five*c2-three)*(five*c2-two)*(five*c2-one)*c2

        shl(9,l)= ((25.d0/12.d0)*(five*c3-one))
     &                     *c3*(five*c2-two)*(five*c2-one)*c2

        shl(10,l)= ((25.d0/12.d0)*(five*c3-two))*(five*c3-one)
     &                      *c3*(five*c2-one)*c2

        shl(11,l)= ((25.d0/24.d0)*(five*c3-three))
     &                     *(five*c3-two)*(five*c3-one)*c3*c2
              
        shl(12,l)= ((25.d0/24.d0)*(five*c3-three))
     &                      *(five*c3-two)*(five*c3-one)*c3*c1
              
        shl(13,l)= ((25.d0/12.d0)*(five*c3-two))
     &                      *(five*c3-one)*c3*(five*c1-one)*c1
              
        shl(14,l)= ((25.d0/12.d0)*(five*c3-one))
     &                       *c3*(five*c1-two)*(five*c1-one)*c1

        shl(15,l)= (25.d0/24.d0)*c3*(five*c1-three)*(five*c1-two)
     &         *(five*c1-one)*c1

        shl(16,l)=((125.d0/6.d0)*(five*c1-two))*(five*c1-one)*c1*c2*c3

        shl(17,l)=((125.d0/4.d0)*(five*c1-one))*(five*c2-one)*c1*c2*c3
              
        shl(18,l)=((125.d0/6.d0)*(five*c2-two))*(five*c2-one)*c1*c2*c3
              
        shl(19,l)= ((125.d0/4.d0)*(five*c2-one))*(five*c3-one)*c1*c2*c3
              
        shl(20,l)= ((125.d0/6.d0)*(five*c3-two))*(five*c3-one)*c1*c2*c3

        shl(21,l)= ((125.d0/4.d0)*(five*c1-one))*(five*c3-one)*c1*c2*c3

            end if
            
  200      continue 
c 

      return
      end
c**** new ********************************************************************** 
      subroutine shapen(shl,nen,nenc)
c
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension raone(8),xaone(8)
      dimension shlone(8,8),inod(8,8),inodc(8,8)
      dimension shl(64,64)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      if (nen.eq.1) xaone(1) = zero
      if (nenc.eq.1) raone(1) = zero
c
c
      if (nen.eq.4) then
         xaone(1) = -one
         xaone(2) =  one
	nenx=2
	neny=2
      endif
c
      if (nenc.eq.4) then
         raone(1) = -one
         raone(2) =  one
	nenrx=2
	nenry=2
      endif
c
c 
      if(nen.eq.9) then
         xaone(1)= -one
         xaone(2)= one
         xaone(3)= zero
	nenx=3
	neny=3
      endif
c 
      if(nenc.eq.9) then
         raone(1)= -one
         raone(2)= one
         raone(3)= zero
	nenrx=3
	nenry=3
      endif
c
c
      if (nen.eq.16) then
         xaone(1) = -one
         xaone(2) = one
         xaone(3) = -.333333333333333
         xaone(4) =  .333333333333333
	nenx=4
	neny=4
         endif
c
      if (nenc.eq.16) then
         raone(1) = -one
         raone(2) = one
         raone(3) = -.333333333333333
         raone(4) =  .333333333333333
	nenrx=4
	nenry=4
         endif
c
c
       if(nen.eq.25) then
         xaone(1)= -one 
         xaone(2)=  one 
         xaone(3)= -pt5
         xaone(4)= zero
         xaone(5)= pt5
	nenx=5
	neny=5
       endif
c
       if(nenc.eq.25) then
         raone(1)= -one 
         raone(2)=  one 
         raone(3)= -pt5
         raone(4)= zero
         raone(5)= pt5
	nenrx=5
	nenry=5
       endif
c
c
        if(nen.eq.36) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.600000000000000
         xaone(4) = -.200000000000000
         xaone(5) =  .200000000000000
         xaone(6) =  .600000000000000
	nenx=6
	neny=6
        endif
c
        if(nenc.eq.36) then
         raone(1) = -one
         raone(2) =  one
         raone(3) = -.600000000000000
         raone(4) = -.200000000000000
         raone(5) =  .200000000000000
         raone(6) =  .600000000000000
	nenrx=6
	nenry=6
        endif
c
c
        if(nen.eq.49) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.666666666666666
         xaone(4) = -.333333333333333
         xaone(5) = zero
         xaone(6) =  .333333333333333
         xaone(7) =  .666666666666666
	nenx=7
	neny=7
        endif
c
        if(nenc.eq.49) then
         raone(1) = -one
         raone(2) =  one
         raone(3) = -.666666666666666
         raone(4) = -.333333333333333
         raone(5) = zero
         raone(6) =  .333333333333333
         raone(7) =  .666666666666666
	nenrx=7
	nenry=7
        endif
c
c
        if(nen.eq.64) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -0.71428571428571
         xaone(4) = -0.42857142857143
         xaone(5) = -0.14285714285714
         xaone(6) = 0.14285714285714
         xaone(7) = 0.42857142857143
         xaone(8) = 0.71428571428571
	nenx=8
	neny=8
        endif
c
        if(nenc.eq.64) then
         raone(1) = -one
         raone(2) =  one
         raone(3) = -0.71428571428571
         raone(4) = -0.42857142857143
         raone(5) = -0.14285714285714
         raone(6) = 0.14285714285714
         raone(7) = 0.42857142857143
         raone(8) = 0.71428571428571
	nenrx=8
	nenry=8
        endif
c
c
	if(nen.eq.1) inod(1,1) = 1
	if(nen.eq.4) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(1,2) = 4
	  inod(2,2) = 3
	end if
c
	if(nen.eq.9) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 7
c
	  inod(1,3) = 8
	  inod(2,3) = 6
	  inod(3,3) = 9
	end if
c
c
c
	if(nen.eq.16) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 10
	  inod(4,2) = 9
c
	  inod(1,3) = 12
	  inod(2,3) = 7
	  inod(3,3) = 13
	  inod(4,3) = 14
c
	  inod(1,4) = 11
	  inod(2,4) = 8
	  inod(3,4) = 16
	  inod(4,4) = 15
	end if
c
c
	if(nen.eq.25) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 13
	  inod(4,2) = 12
	  inod(5,2) = 11
c
	  inod(1,3) = 16
	  inod(2,3) = 8
	  inod(3,3) = 17
	  inod(4,3) = 18
	  inod(5,3) = 19
c
	  inod(1,4) = 15
	  inod(2,4) = 9
	  inod(3,4) = 24
	  inod(4,4) = 25
	  inod(5,4) = 20
c	  
	  inod(1,5) = 14
	  inod(2,5) = 10
	  inod(3,5) = 23
	  inod(4,5) = 22
	  inod(5,5) = 21
	end if
c
c
	if(nen.eq.36) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 16
	  inod(4,2) = 15
	  inod(5,2) = 14
	  inod(6,2) = 13
c
	  inod(1,3) = 20
	  inod(2,3) = 9
	  inod(3,3) = 21
	  inod(4,3) = 22
	  inod(5,3) = 23
	  inod(6,3) = 24
c
	  inod(1,4) = 19
	  inod(2,4) = 10
	  inod(3,4) = 32
	  inod(4,4) = 33
	  inod(5,4) = 34
	  inod(6,4) = 25
c	  
	  inod(1,5) = 18
	  inod(2,5) = 11
	  inod(3,5) = 31
	  inod(4,5) = 36
	  inod(5,5) = 35
	  inod(6,5) = 26
c
	  inod(1,6) = 17
	  inod(2,6) = 12
	  inod(3,6) = 30
	  inod(4,6) = 29
	  inod(5,6) = 28
	  inod(6,6) = 27
c	  
	end if
c
c
	if(nen.eq.49) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
	  inod(7,1) = 9
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 19
	  inod(4,2) = 18
	  inod(5,2) = 17
	  inod(6,2) = 16
	  inod(7,2) = 15
c
	  inod(1,3) = 24
	  inod(2,3) = 10
	  inod(3,3) = 25
	  inod(4,3) = 26
	  inod(5,3) = 27
	  inod(6,3) = 28
	  inod(7,3) = 29
c
	  inod(1,4) = 23
	  inod(2,4) = 11
	  inod(3,4) = 40
	  inod(4,4) = 41
	  inod(5,4) = 42
	  inod(6,4) = 43
	  inod(7,4) = 30
c	  
	  inod(1,5) = 22
	  inod(2,5) = 12
	  inod(3,5) = 39
	  inod(4,5) = 48
	  inod(5,5) = 49
	  inod(6,5) = 44
	  inod(7,5) = 31
c
	  inod(1,6) = 21
	  inod(2,6) = 13
	  inod(3,6) = 38
	  inod(4,6) = 47
	  inod(5,6) = 46
	  inod(6,6) = 45
	  inod(7,6) = 32
c
	  inod(1,7) = 20
	  inod(2,7) = 14
	  inod(3,7) = 37
	  inod(4,7) = 36
	  inod(5,7) = 35
	  inod(6,7) = 34
	  inod(7,7) = 33
c	  
	end if
c
c
	if(nen.eq.64) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
	  inod(7,1) = 9
	  inod(8,1) = 10
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 22
	  inod(4,2) = 21
	  inod(5,2) = 20
	  inod(6,2) = 19
	  inod(7,2) = 18
	  inod(8,2) = 17
c
	  inod(1,3) = 28
	  inod(2,3) = 11
	  inod(3,3) = 29
	  inod(4,3) = 30
	  inod(5,3) = 31
	  inod(6,3) = 32
	  inod(7,3) = 33
	  inod(8,3) = 34
c
	  inod(1,4) = 27
	  inod(2,4) = 12
	  inod(3,4) = 48
	  inod(4,4) = 49
	  inod(5,4) = 50
	  inod(6,4) = 51
	  inod(7,4) = 52
	  inod(8,4) = 35
c	  
	  inod(1,5) = 26
	  inod(2,5) = 13
	  inod(3,5) = 47
	  inod(4,5) = 60
	  inod(5,5) = 61
	  inod(6,5) = 62
	  inod(7,5) = 53
	  inod(8,5) = 36
c
	  inod(1,6) = 25
	  inod(2,6) = 14
	  inod(3,6) = 46
	  inod(4,6) = 59
	  inod(5,6) = 64
	  inod(6,6) = 63
	  inod(7,6) = 54
	  inod(8,6) = 37
c
	  inod(1,7) = 24
	  inod(2,7) = 15
	  inod(3,7) = 45
	  inod(4,7) = 58
	  inod(5,7) = 57
	  inod(6,7) = 56
	  inod(7,7) = 55
	  inod(8,7) = 38
c
	  inod(1,8) = 23
	  inod(2,8) = 16
	  inod(3,8) = 44
	  inod(4,8) = 43
	  inod(5,8) = 42
	  inod(6,8) = 41
	  inod(7,8) = 40
	  inod(8,8) = 39
c	  
	end if
c
c
c
	if(nenc.eq.1) inodc(1,1) = 1
	if(nenc.eq.4) then
	  inodc(1,1) = 1
	  inodc(2,1) = 2
	  inodc(1,2) = 4
	  inodc(2,2) = 3
	end if
c
	if(nenc.eq.9) then
	  inodc(1,1) = 1
	  inodc(2,1) = 2
	  inodc(3,1) = 5
c
	  inodc(1,2) = 4
	  inodc(2,2) = 3
	  inodc(3,2) = 7
c
	  inodc(1,3) = 8
	  inodc(2,3) = 6
	  inodc(3,3) = 9
	end if
c
c
c
	if(nenc.eq.16) then
	  inodc(1,1) = 1
	  inodc(2,1) = 2
	  inodc(3,1) = 5
	  inodc(4,1) = 6
c
	  inodc(1,2) = 4
	  inodc(2,2) = 3
	  inodc(3,2) = 10
	  inodc(4,2) = 9
c
	  inodc(1,3) = 12
	  inodc(2,3) = 7
	  inodc(3,3) = 13
	  inodc(4,3) = 14
c
	  inodc(1,4) = 11
	  inodc(2,4) = 8
	  inodc(3,4) = 16
	  inodc(4,4) = 15
	end if
c
c
	if(nenc.eq.25) then
	  inodc(1,1) = 1
	  inodc(2,1) = 2
	  inodc(3,1) = 5
	  inodc(4,1) = 6
	  inodc(5,1) = 7
c
	  inodc(1,2) = 4
	  inodc(2,2) = 3
	  inodc(3,2) = 13
	  inodc(4,2) = 12
	  inodc(5,2) = 11
c
	  inodc(1,3) = 16
	  inodc(2,3) = 8
	  inodc(3,3) = 17
	  inodc(4,3) = 18
	  inodc(5,3) = 19
c
	  inodc(1,4) = 15
	  inodc(2,4) = 9
	  inodc(3,4) = 24
	  inodc(4,4) = 25
	  inodc(5,4) = 20
c	  
	  inodc(1,5) = 14
	  inodc(2,5) = 10
	  inodc(3,5) = 23
	  inodc(4,5) = 22
	  inodc(5,5) = 21
	end if
c
c
	if(nenc.eq.36) then
	  inodc(1,1) = 1
	  inodc(2,1) = 2
	  inodc(3,1) = 5
	  inodc(4,1) = 6
	  inodc(5,1) = 7
	  inodc(6,1) = 8
c
	  inodc(1,2) = 4
	  inodc(2,2) = 3
	  inodc(3,2) = 16
	  inodc(4,2) = 15
	  inodc(5,2) = 14
	  inodc(6,2) = 13
c
	  inodc(1,3) = 20
	  inodc(2,3) = 9
	  inodc(3,3) = 21
	  inodc(4,3) = 22
	  inodc(5,3) = 23
	  inodc(6,3) = 24
c
	  inodc(1,4) = 19
	  inodc(2,4) = 10
	  inodc(3,4) = 32
	  inodc(4,4) = 33
	  inodc(5,4) = 34
	  inodc(6,4) = 25
c	  
	  inodc(1,5) = 18
	  inodc(2,5) = 11
	  inodc(3,5) = 31
	  inodc(4,5) = 36
	  inodc(5,5) = 35
	  inodc(6,5) = 26
c
	  inodc(1,6) = 17
	  inodc(2,6) = 12
	  inodc(3,6) = 30
	  inodc(4,6) = 29
	  inodc(5,6) = 28
	  inodc(6,6) = 27
c	  
	end if
c
c
	if(nenc.eq.49) then
	  inodc(1,1) = 1
	  inodc(2,1) = 2
	  inodc(3,1) = 5
	  inodc(4,1) = 6
	  inodc(5,1) = 7
	  inodc(6,1) = 8
	  inodc(7,1) = 9
c
	  inodc(1,2) = 4
	  inodc(2,2) = 3
	  inodc(3,2) = 19
	  inodc(4,2) = 18
	  inodc(5,2) = 17
	  inodc(6,2) = 16
	  inodc(7,2) = 15
c
	  inodc(1,3) = 24
	  inodc(2,3) = 10
	  inodc(3,3) = 25
	  inodc(4,3) = 26
	  inodc(5,3) = 27
	  inodc(6,3) = 28
	  inodc(7,3) = 29
c
	  inodc(1,4) = 23
	  inodc(2,4) = 11
	  inodc(3,4) = 40
	  inodc(4,4) = 41
	  inodc(5,4) = 42
	  inodc(6,4) = 43
	  inodc(7,4) = 30
c	  
	  inodc(1,5) = 22
	  inodc(2,5) = 12
	  inodc(3,5) = 39
	  inodc(4,5) = 48
	  inodc(5,5) = 49
	  inodc(6,5) = 44
	  inodc(7,5) = 31
c
	  inodc(1,6) = 21
	  inodc(2,6) = 13
	  inodc(3,6) = 38
	  inodc(4,6) = 47
	  inodc(5,6) = 46
	  inodc(6,6) = 45
	  inodc(7,6) = 32
c
	  inodc(1,7) = 20
	  inodc(2,7) = 14
	  inodc(3,7) = 37
	  inodc(4,7) = 36
	  inodc(5,7) = 35
	  inodc(6,7) = 34
	  inodc(7,7) = 33
c	  
	end if
c
c
	if(nenc.eq.64) then
	  inodc(1,1) = 1
	  inodc(2,1) = 2
	  inodc(3,1) = 5
	  inodc(4,1) = 6
	  inodc(5,1) = 7
	  inodc(6,1) = 8
	  inodc(7,1) = 9
	  inodc(8,1) = 10
c
	  inodc(1,2) = 4
	  inodc(2,2) = 3
	  inodc(3,2) = 22
	  inodc(4,2) = 21
	  inodc(5,2) = 20
	  inodc(6,2) = 19
	  inodc(7,2) = 18
	  inodc(8,2) = 17
c
	  inodc(1,3) = 28
	  inodc(2,3) = 11
	  inodc(3,3) = 29
	  inodc(4,3) = 30
	  inodc(5,3) = 31
	  inodc(6,3) = 32
	  inodc(7,3) = 33
	  inodc(8,3) = 34
c
	  inodc(1,4) = 27
	  inodc(2,4) = 12
	  inodc(3,4) = 48
	  inodc(4,4) = 49
	  inodc(5,4) = 50
	  inodc(6,4) = 51
	  inodc(7,4) = 52
	  inodc(8,4) = 35
c	  
	  inodc(1,5) = 26
	  inodc(2,5) = 13
	  inodc(3,5) = 47
	  inodc(4,5) = 60
	  inodc(5,5) = 61
	  inodc(6,5) = 62
	  inodc(7,5) = 53
	  inodc(8,5) = 36
c
	  inodc(1,6) = 25
	  inodc(2,6) = 14
	  inodc(3,6) = 46
	  inodc(4,6) = 59
	  inodc(5,6) = 64
	  inodc(6,6) = 63
	  inodc(7,6) = 54
	  inodc(8,6) = 37
c
	  inodc(1,7) = 24
	  inodc(2,7) = 15
	  inodc(3,7) = 45
	  inodc(4,7) = 58
	  inodc(5,7) = 57
	  inodc(6,7) = 56
	  inodc(7,7) = 55
	  inodc(8,7) = 38
c
	  inodc(1,8) = 23
	  inodc(2,8) = 16
	  inodc(3,8) = 44
	  inodc(4,8) = 43
	  inodc(5,8) = 42
	  inodc(6,8) = 41
	  inodc(7,8) = 40
	  inodc(8,8) = 39
c	  
	end if
c
      do 100 l = 1, nenrx
         r = raone(l)
c
        if(nenx.eq.1) then
         shlone(1,l) = zero
         go to 100
        endif
c
        do 50 i = 1, nenx
         aa = one
         bb = one
         aax = zero
         do 40 j =1, nenx
          daj = one
          if (i .ne. j)then
          aa = aa * ( r - xaone(j))
          bb = bb * ( xaone(i) - xaone(j))
          endif
   40    continue
        shlone(i,l) = aa/bb
   50  continue
c
  100  continue      
c
      do ly=1,nenrx
      do lx=1,nenry
        l = inodc(lx,ly)
	 do iy=1,neny
	 do ix=1,nenx
         j = inod(ix,iy)
	   shl(j,l) = shlone(ix,lx)*shlone(iy,ly)
	 end do
	 end do
	end do
	end do
c
      return
      end
      
c***********************************************************************
c          resolução do sistema
c***********************************************************************
c**** new **********************************************************************
      subroutine formlm (id,ien,lm,ndof,ned,nen,numel)
c
c.... program to form lm array
c
      dimension id(ndof,*),ien(nen,*),lm(ned,nen,*)
c
      do 300 k=1,numel
c
      do 200 j=1,nen
      node=ien(j,k)
c
      do 100 i=1,ndof
      lm(i,j,k) = id(i,node)
  100 continue
c
  200 continue
c
  300 continue
c
      return
      end
c**** new **********************************************************************
      subroutine colht(idiag,lm,ned,nen,numel,neq)
c
c.... program to compute column heights in global left-hand-side matrix
c
      dimension idiag(*),lm(ned,nen,*)
c
      do 500 k=1,numel
      min = neq
c
      do 200 j=1,nen
c
      do 100 i=1,ned
      num = lm(i,j,k)
      if (num.gt.0) min = min0(min,num)
  100 continue
c
  200 continue
c
      do 400 j=1,nen
c
      do 300 i=1,ned
      num = lm(i,j,k)
      if (num.gt.0) then
         m = num - min
         if (m.gt.idiag(num)) idiag(num) = m
      endif
c
  300 continue

  400 continue
c
  500 continue
c
      return
      end
c**** new **********************************************************************
      subroutine diag(idiag,neq,n)
c
c.... program to compute diagonal addresses of left-hand-side matrix
c
      dimension idiag(*)
c
      n = 1
      idiag(1) = 1 
      if (neq.eq.1) return
c
      do 100 i=2,neq
      idiag(i) = idiag(i) + idiag(i-1) + 1
  100 continue
      n = idiag(neq)
c
      return
      end
      
c************************************************************************
       subroutine condtdg(elmbb,elmcb,elmdb,
     &            elfbb,eleffd,elresd,
     &            neep,nee)
c
      implicit real*8 (a-h,o-z)
c                                                                       
      dimension elmbb(neep,*),elmcb(neep,*),elmdb(nee,*),
     &          eleffd(nee,*)
      dimension elfbb(*),elresd(*)
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
c
c   condensacao do sitema
c
c
c   Bb Xb   + Cb Xc    = Fbb
c
c   Cb^t Xb  + Cc Xc   = Fcb
c
c   com eliminacao das incognitas  Xb
c
      call invmb(elmbb,neep,neep)            
c
c    \br(D) = Cb^t Bb^{1}
c
      do i=1,nee
      do j=1,neep
        elmdb(i,j) = 0.d00
        do k=1,neep
          elmdb(i,j) = elmdb(i,j) + elmcb(k,i)*elmbb(k,j)
        end do
      end do
      end do
c
c  eleffd = eleffd - C^t B^{1}C
c
      do i=1,nee
      do j=1,nee
        do k=1,neep
          eleffd(i,j) = eleffd(i,j) - elmdb(i,k)*elmcb(k,j)
        end do
      end do
      end do
c
c  elresd = elresd - C^t B^{1}Fb
c
      do i=1,nee
        do k=1,neep
          elresd(i) = elresd(i) - elmdb(i,k)*elfbb(k)
        end do
      end do

      return
c
      end
c**** new **********************************************************************
       subroutine condtdg2(elmcc,eleffd,elfcc,elresd,nee)
c
      implicit real*8 (a-h,o-z)
c
      dimension elmcc(nee+1,*),elfcc(*)
      dimension eleffd(nee,*),elresd(*)
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
c
c   A Xu    +  B Xp  = Fb
c
c   B^t Xu  +  C Xp  = Fc
c
c   com eliminacao das incognitas  Xp
c    ______________
c   |           |  |
c   |           |  |
c   |           |  |
c   |     A     |B |
c   |           |  |      Matriz eleffd
c   |           |  |
c   |___________|__|
c   |    B^t    |C |
c   |___________|__|
c
c
c   Xp = r*(Fc - B^t Xu)
c
c
c ....r = C^{-1} 
c
      r = 1.d00/(elmcc(nee+1,nee+1)) 
c 
c ....A = A - r * B * B^t    
c
      do i=1,nee 
         elresd(i) = elfcc(i)
      do j=1,nee
	   eleffd(i,j) = elmcc(i,j) 
	end do 
	end do
      do i=1,nee 
      do j=1,nee
	   eleffd(i,j) = eleffd(i,j) 
     &             - r*elmcc(i,nee+1)*elmcc(nee+1,j)
	  end do 
	  end do
c
c     Fb = Fb - r * B * Fc
c
      do i=1,nee
       elresd(i) = elresd(i) 
     &      	 - r*elmcc(i,nee+1)*elfcc(nee+1)
	  end do 
c

      return
c
      end
c**** new ********************************************************************** 
       subroutine condntdg(elmbb,elmcb,elmbc,elmdb,
     &            elfbb,eleffd,elresd,
     &            neep,nee)
c
      implicit real*8 (a-h,o-z)
c                                                                       
      dimension elmbb(neep,*),elmcb(neep,*),elmdb(nee,*),
     &          eleffd(nee,*),elmbc(nee,*)
      dimension elfbb(*),elresd(*)
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2,ierrb0,ierrb,
     & ierrbi
c
c
c   condensa??o do sitema
c
c
c   Bb Xb   + Cb Xc    = Fbb
c
c   Bc Xb  + Cc Xc   = Fc
c
c   com elimina??o das inciggnitas  Xb
c
      call invmb(elmbb,neep,neep)            
c
c    \br(D) = Bc * B^{1}
c
      do i=1,nee
      do j=1,neep
        elmdb(i,j) = 0.d00
        do k=1,neep
          elmdb(i,j) = elmdb(i,j) + elmbc(i,k)*elmbb(k,j)
        end do
      end do
      end do
c
c  eleffd = eleffd - Bc B^{1}C
c
      do i=1,nee
      do j=1,nee
        do k=1,neep
          eleffd(i,j) = eleffd(i,j) - elmdb(i,k)*elmcb(k,j)
        end do
      end do
      end do
c
c  elresd = elresd - Bc B^{1}Fb
c
      do i=1,nee
        do k=1,neep
          elresd(i) = elresd(i) - elmdb(i,k)*elfbb(k)
        end do
      end do

      return
c
      end

c**** new **********************************************************************
      subroutine invmb(am,ndim,m)                      
c                                                                        
c     subrotina de inversao                                              
c                                                                        
      implicit real*8(a-h,o-z)                                  
      dimension dis(200,1),piv(200),am(ndim,*),ipi(200),ind(200,2)
      ncoln=0                                                            
      det=1.0d00                                                            
      do 20 j=1,m                                                        
   20 ipi(j)=0                                                           
      do 550 i=1,m                                                       
      amax=0.0d00                                                           
      do 105 j=1,m                                                       
      if(ipi(j)-1)60,105,60                                              
   60 do 100 k=1,m                                                       
      if(ipi(k)-1) 80,100,740                                            
   80 if(dabs(amax)-dabs(am(j,k)))85,100,100                             
   85 irow=j                                                             
      ico=k                                                              
      amax=am(j,k)                                                       
  100 continue                                                           
  105 continue                                                           
      ipi(ico)=ipi(ico)+1                                                
      if(irow-ico)140,260,140                                            
  140 det=-det                                                           
      do 200 l=1,m                                                       
      swap=am(irow,l)                                                    
      am(irow,l)=am(ico,l)                                               
  200 am(ico,l)=swap                                                     
      if(ncoln) 260,260,210                                              
  210 do 250 l=1,ncoln                                                   
      swap=dis(irow,l)                                                   
      dis(irow,l)=dis(ico,l)                                             
  250 dis(ico,l)=swap                                                    
  260 ind(i,1)=irow                                                      
      piv(i)=am(ico,ico)                                                 
      ind(i,2)=ico                                                       
      det=det*piv(i)                                                     
      am(ico,ico)=1.0d00                                                    
      do 350 l=1,m                                                       
  350 am(ico,l)=am(ico,l)/piv(i)                                         
      if(ncoln) 380,380,360                                              
  360 do 370 l=1,ncoln                                                   
  370 dis(ico,l)=dis(ico,l)/piv(i)                                       
  380 do 550 lz=1,m                                                      
      if(lz-ico)400,550,400                                              
  400 t=am(lz,ico)                                                       
      am(lz,ico)=0.0d00                                                     
      do 450 l=1,m                                                       
  450 am(lz,l)=am(lz,l)-am(ico,l)*t                                      
      if(ncoln)550,550,460                                               
  460 do 500 l=1,ncoln                                                   
  500 dis(lz,l)=dis(lz,l)-dis(ico,l)*t                                   
  550 continue                                                           
      do 710 i=1,m                                                       
      l=m+1-i                                                            
      if(ind(l,1)-ind(l,2))630,710,630                                   
  630 jrow=ind(l,1)                                                      
      jco=ind(l,2)                                                       
      do 705 k=1,m                                                       
      swap=am(k,jrow)                                                    
      am(k,jrow)=am(k,jco)                                               
      am(k,jco)=swap                                                     
  705 continue                                                           
  710 continue                                                           
  740 continue                                                           
      return                                                             
      end                                                                
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine addnsl(alhs,clhs,eleffm,idiag,lm,nee,ldiag)
c
c         program to add element left-hand-side matrix to          
c                global left-hand-side matrix                      
c                                                                  
c        ldiag = .true.,  add diagonal element matrix              
c                                                                  
c        ldiag = .false, then                                     
c        add full nonsymmetric element matrix                   
c                                                                  
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      logical ldiag
      dimension alhs(*),clhs(*),eleffm(nee,*),idiag(*),lm(*)
c
      if (ldiag) then
c
         do 100 j=1,nee
            k = iabs(lm(j))
            if (k.gt.0) then
               l = idiag(k)
               alhs(l) = alhs(l) + eleffm(j,j)
            endif
  100    continue
c
      else
c
         do 400 j=1,nee
            k = iabs(lm(j))
            if (k.gt.0) then
               do 200 i=1,nee
                  m = iabs(lm(i))
                  if (m.gt.0) then
                     if (k.gt.m) then
                        l = idiag(k) - k + m
                        alhs(l) = alhs(l) + eleffm(i,j)
                     else
                        l = idiag(m) - m + k
                        clhs(l) = clhs(l) + eleffm(i,j)
                     endif
                     if (k.eq.m) then
                        l = idiag(k)
                        alhs(l) = alhs(l) + eleffm(i,j)
                        clhs(l) = alhs(l)
                     endif
                  endif
  200          continue
            endif
  400    continue
c
      endif
c
      return
      end

c**** new **********************************************************************
      subroutine addlhs(alhs,eleffm,idiag,lm,nee,diag)
c
c.... program to add element left-hand-side matrix to
c        global left-hand-side matrix
c
c        diag = .true., add diagonal element matrix
c
c        diag = .false, add upper triangle of full element matrix
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      logical diag
      dimension alhs(*),eleffm(nee,*),idiag(*),lm(*)
c
      if (diag) then
c
         do 100 j=1,nee
         k = lm(j)
         if (k.gt.0) then
            l = idiag(k)
            alhs(l) = alhs(l) + eleffm(j,j)
         endif
  100    continue
c
      else
c
         do 300 j=1,nee
         k = lm(j)
         if (k.gt.0) then
c
            do 200 i=1,j
            m = lm(i)
            if (m.gt.0) then
               if (k.ge.m) then
                  l = idiag(k) - k + m
               else
                  l = idiag(m) - m + k
               endif
               alhs(l) = alhs(l) + eleffm(i,j)
            endif
  200       continue
c
         endif
  300    continue
c
      endif
c
      return
      end
c**** new **********************************************************************
      subroutine addrhs(brhs,elresf,lm,nee)
c
c.... program to add element residual-force vector to
c        global right-hand-side vector
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension brhs(*),elresf(*),lm(*)
c
      do 100 j=1,nee
      k = lm(j)
      if (k.gt.0) brhs(k) = brhs(k) + elresf(j)
  100 continue
c
      return
      end
c**** newfrey *********************************************************
      subroutine factns(a,c,idiag,neq)
c
c.... program to perform crout factorization: a = l * d * u
c
c        a(i):  coefficient matrix stored in compacted column form;
c               after factorization contains d and u
c
c        c(i):  non-symmetric lower triangular coefficient matrix stored in
c                compacted row form; after factorization contains l
c 
c
      implicit real*8 (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(*),c(*),idiag(*)
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      jj = 0
c
      do 300 j=1,neq
c
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
c
      if (jcolht.gt.2) then
c
c....... for column j and i.le.j-1, replace a(i,j) with d(i,i)*u(i,j)
c
         istart = j - jcolht + 2
         jm1    = j - 1
         ij     = jjlast + 2
         ii     = idiag(istart-1)
c
         do 100 i=istart,jm1
c
         iilast = ii
         ii     = idiag(i)
         icolht = ii - iilast
         length = min0(icolht-1,i - istart + 1)
         if (length.gt.0)  then
            a(ij) = a(ij) - coldot(a(ij-length),c(ii-length),length)
            c(ij) = c(ij) - coldot(c(ij-length),a(ii-length),length)
         endif
         ij = ij + 1
  100    continue
c
      endif
c
      if (jcolht.ge.2) then
c
c....... for column j and i.le.j-1, replace a(i,j) with u(i,j);
c           replace a(j,j) with d(j,j).
c
         jtemp = j - jj
c
         do 200 ij=jjlast+1,jj-1
c
         ii = idiag(jtemp + ij)
c
c....... warning: the following calculations are skipped 
c                 if a(ii) equals zero
c
         if (a(ii).ne.zero) then
             c(ij) = c(ij)/a(ii)
             a(jj) = a(jj) - c(ij)*a(ij)
             a(ij) = a(ij)/a(ii)
         endif
  200    continue
c
      endif
c
  300 continue
c
      return
      end
c**** new **********************************************************************
      subroutine factor(a,idiag,neq)
c
c.... program to perform Crout factorization: a = u(transpose) * d * u
c
c        a(i):  coefficient matrix stored in compacted column form;
c               after factorization contains d and u
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension a(*),idiag(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      jj = 0
c
      do 300 j=1,neq
c
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
c
      if (jcolht.gt.2) then
c
c....... for column j and i.le.j-1, replace a(i,j) with d(i,i)*u(i,j)
c
         istart = j - jcolht + 2
         jm1    = j - 1
         ij     = jjlast + 2
         ii     = idiag(istart-1)
c
         do 100 i=istart,jm1
c
         iilast = ii
         ii     = idiag(i)
         icolht = ii - iilast
         jlngth = i - istart + 1
         length = min0(icolht-1,jlngth)
         if (length.gt.0) 
     &      a(ij) = a(ij) - coldot(a(ii-length),a(ij-length),length)
         ij = ij + 1
  100    continue
c
      endif
c
      if (jcolht.ge.2) then
c
c....... for column j and i.le.j-1, replace a(i,j) with u(i,j);
c           replace a(j,j) with d(j,j).
c
         jtemp = j - jj
c
         do 200 ij=jjlast+1,jj-1
c
         ii = idiag(jtemp + ij)
         if (a(ii).ne.zero) then
            temp  = a(ij)
            a(ij) = temp/a(ii)
            a(jj) = a(jj) - temp*a(ij)
         endif
  200    continue
c
      endif
c
  300 continue
c
      return
      end
      
c**** new **********************************************************************
      subroutine solvetdg(elmbb,elfbb,neep)
c
      implicit real*8 (a-h,o-z)
c                                                                       
      dimension elmbb(neep,*),elfbb(*)
      dimension fab(200)
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
c   resolve o sistema
c
c   B Xb = Fb
c
      call invmb(elmbb,neep,neep)     
c
c    Aa = A^{-1} Fa
c 
      do j=1,neep
        fab(j)=0.d00
        do k=1,neep
          fab(j) = fab(j) + elmbb(j,k)*elfbb(k)
        end do
      end do
c
c  xb=fab
c
      do i=1,neep
      elfbb(i) = fab(i)
      end do      
c
      return
c
      end
c**** new **********************************************************************
      subroutine back(a,b,idiag,neq)
c
c.... program to perform forward reduction and back substitution
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension a(*),b(*),idiag(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
c.... forward reduction
c
      jj = 0
c
      do 100 j=1,neq
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
      if (jcolht.gt.1) 
     &   b(j) = b(j) - coldot(a(jjlast+1),b(j-jcolht+1),jcolht-1)
  100 continue
c
c.... diagonal scaling
c
      do 200 j=1,neq
      ajj = a(idiag(j))
      if (ajj.ne.zero) b(j) = b(j)/ajj
  200 continue
c
c.... back substitution
c
      if (neq.eq.1) return
      jjnext = idiag(neq)
c
      do 400 j=neq,2,-1
      jj     = jjnext
      jjnext = idiag(j-1)
      jcolht = jj - jjnext
      if (jcolht.gt.1) then
         bj = b(j)
         istart = j - jcolht + 1
         jtemp  = jjnext - istart + 1
c
         do 300 i=istart,j-1
         b(i) = b(i) - a(jtemp+i)*bj
  300    continue
c
      endif
c
  400 continue
c
      return
      end
c**** newfrey *********************************************************
      subroutine backns(a,c,b,idiag,neq)
c
c.... program to perform forward reduction and back substitution
c
      implicit real*8 (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(*),c(*),b(*),idiag(*)
c
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
c.... forward reduction
c
c
      jj = 0
c
      do 100 j=1,neq
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
      if (jcolht.gt.1) then
          b(j) = b(j) - coldot(c(jjlast+1),b(j-jcolht+1),jcolht-1)
      endif
  100 continue
c
c.... diagonal scaling
c
      do 200 j=1,neq
      ajj = a(idiag(j))
c
c.... warning: diagonal scaling is not performed if ajj equals zero
c
      if (ajj.ne.zero) b(j) = b(j)/ajj
  200 continue
c
c.... back substitution
c
      if (neq.eq.1) return
      jjnext = idiag(neq)
c
      do 400 j=neq,2,-1
      jj     = jjnext
      jjnext = idiag(j-1)
      jcolht = jj - jjnext
      if (jcolht.gt.1) then
         bj = b(j)
         istart = j - jcolht + 1
         jtemp  = jjnext - istart + 1
c
         do 300 i=istart,j-1
         b(i) = b(i) - a(jtemp+i)*bj
  300    continue
c
      endif
c
  400 continue
c
      return
      end
      
c**** new **********************************************************************
      subroutine btod(id,d,brhs,ndof,numnp)
c
c.... program to perform transfer from r.h.s. to displacement array
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension id(ndof,*),d(ndof,*),brhs(*)
c
         do 200 i=1,ndof
c
         do 100 j=1,numnp
         k = id(i,j)
         if (k.gt.0) d(i,j) = brhs(k)
  100    continue
c
  200    continue
c
      return
      end
c**** new **********************************************************************
      subroutine btodlocal(d,brhs,ndof,numnp,nel)
c
c.... program to perform transfer from local r.h.s. to local displacement array
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension d(ndof,numnp,*),brhs(*)
c
      do  j=1,numnp 
         do  i=1,ndof
            k=(j-1)*ndof + i
            d(i,j,nel) = brhs(k)
         end do
      end do 
c
      return
      end
c***********************************************************************
c          estimativas de erro
c***********************************************************************
      
c**** new **********************************************************************
      subroutine flnorm(ien   ,x     ,xl   ,    
     &                 d     ,dl    ,mat   ,
     &                 matside,
     &                 c     ,ipar  ,dlf   ,   
     &                 dlp   ,dsfl  ,det   ,
     &                 shl   , shg  ,wt    ,
     &                 detc  ,shlc  , shgc , 
     &                 ddis  ,detp  ,shlp  , 
     &                 shgp  ,lado  ,
c     
     &                 ideg  ,
c
     &                 shln  ,shgn  ,
     &                 detn  ,shlb  ,shgb  ,
     &                 detpn ,shlpn ,shgpn ,
     &                 idside,xls   ,idlsd ,
     &                 grav  ,wn    ,
c                  
     &                 numel ,neesq ,nen   ,nsd   ,
     &                 nesd  ,nint  ,neg   ,nrowsh,
     &                 ned   ,nee   ,numnp ,ndof  ,
     &                 ncon  ,nencon,necon ,index ,
     &                 nints ,iwrite ,iplt ,ierb  ,
     &                 nenp  ,nside,nnods ,nenlad, 
     &                 npars ,nmultp, nodsp)
c------------------------------------------------------------------
c     Some variables used in this subroutine
c            j: degree of freedom (1,...ncon)
c         u(j): finite element solution
c        du(j): derivative of finite element solution
c        ue(j): exact solution
c       due(j): derivative of exact solution
c       el2(j): error in L2
c      epri(j): error in the seminorm of H1 (L2 of derivatives)
c     el2el(j): error in L2 in the element domain
c    epriel(j): error in the seminorm of H1 in the element domain
c----------------------------------------------------------------
c
c
c    Program to calculate and print the L2 and H1 seminorm of the error
c    for each degree of freedom in the finite element solution .
c    The trapezoidal rule is used for integration in each element .
c    The number of integration points is given by nints .
c
c    This version is only applicable to 2D.
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      character *1 tab
      dimension ien(nen,*),x(nsd,*),xl(nesd,*),d(ndof,*),dl(ndof,*),
     &          mat(*),c(9,*),ipar(nodsp,*),dlf(ncon,*),matside(*),
     &          dsfl(ncon,nencon,*),ddis(ned,nenp,*)
      dimension u(6),dux(6),duy(6),ue(6),duex(6),duey(6)
      dimension shl(3,nen,*),shg(3,nen,*),det(*),wt(*)
      dimension shlc(3,nencon,*),shgc(3,nencon,*),detc(*)
      dimension el2(6),epri(6),eprix(6),epriy(6),el2el(6),
     &          eprxel(6),epryel(6)
	  dimension detp(*),shlp(3,nenp,*),shgp(3,nenp,*),dlp(ned,*)
c
      dimension detn(*),shlb(3,nenlad,*),shgb(3,nenlad,*),
     &          detpn(*),shlpn(3,npars,*),shgpn(3,npars,*),
     &          idside(nside,*),xls(nesd,*),idlsd(*)
	  dimension dls(ndof,16),grav(*),wn(*)
      dimension shln(3,nnods,*),shgn(3,nnods,*)
c
      dimension ideg(2*ndof,*),lado(nside,*)
c
      common /colhtc/ neq
      common /consts/zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
      pi=4.*datan(1.d00)
      dpi=2.d00*pi
c
      gf1=grav(1)
      gf2=grav(2)
      gf3=grav(3)
      gf4=grav(4)
      gf5=grav(5) 
      gf6=grav(6) 
	  tab=char(9)
c
      xint=0.d00 
      yint=0.d00
c
c    erro na norma l2 da aproximacao, do gradiente e do divergente da velocidade
c
	  edu = 0.d00
	  edux= 0.d00
      eduy= 0.d00
	  divu= 0.d00        
c
c    erro na norma l2 da aproximacao e do gradiente da pressao
c        
	  edp = 0.d00
	  edpx= 0.d00
	  edpy= 0.d00
c
c
c    erro na norma l2 dos multiplicadores (velocidade e pressao)
c
	  xmlt= 0.d00 ! multiplicador da velocidade
c
      call aproxprn(dsfl, ddis, x, ien, mat,
     1                    ncon, nencon, numel, nenp,
     2                    nsd, nen, numnp, neg, numat,
     3                    nesd, iecho,c)
c
c
c.....  loop on elements
c      
      do 50 n=1,numel
c
      mm = mat(n)
        alpha=c(9,mm)
	  xka1=c(1,mm)*alpha
	  xka2=c(2,mm)*alpha
	  eps=c(3,mm)  !
c
                xkl = c(1,mm)
c      
      call local(ien(1,n),x,xl,nen,nsd,nesd)
      call local(ipar(1,n),d,dl,nodsp,ndof,ndof)
c
      do i=1,nenp
	   dlp(1,i) = ddis(1,i,n)
	  end do
c
c
      do 9999 i=1,nencon
        do 8888 j=1,ncon
          dlf(j,i)=dsfl(j,i,n)
8888    continue
9999  continue
c
          un1  = 0.d00
          un2  = 0.d00
          du1x = 0.d00
		  du2x = 0.d00
          du1y = 0.d00
		  du2y = 0.d00
c
		  divue=0.d00
		  dpe  = 0.d00
	      dpex = 0.d00
		  dpey = 0.d00
c
c.....  loop on integration points
c
c
c.....Triangles or quadrilaterals
c
c
ccx       call shgq(xl,det,shl,shg,nint,n,neg,.true.,nen)
c
      call shgqs(xl,detc,shlc,shgc,nint,n,neg,.true.,nencon,shl,nen)
c
      call shgqs(xl,detp,shlp,shgp,nint,n,neg,.true.,nenp,shl,nen)
c
        do 4040 l=1,nint
        ct=detp(l)*wt(l)
		call clear ( u, ncon )
        call clear (dux, ncon )
        call clear (duy, ncon )
        xint=0.d0
        yint=0.d0
         do 3030 i=1,nen
            xint=xint+shl(3,i,l)*xl(1,i)
            yint=yint+shl(3,i,l)*xl(2,i)
3030     continue
c
          do 3330 i=1,nencon 
            do 2020 j=1,ncon
                u(j)   = u(j) + shgc(3,i,l)*dlf(j,i)
	          dux(j) = dux(j) + shgc(1,i,l)*dlf(j,i)
	          duy(j) = duy(j) + shgc(2,i,l)*dlf(j,i)
2020        continue
3330      continue
c
         call uexafx(xint,yint,ue,duex,duey,alpha,mm,index)
c
c    velocidade   
c
         un1  = un1  + ct * ( (u(1)-ue(1))**2 )
         un2  = un2  + ct * ( (u(2)-ue(2))**2 )
         du1x = du1x + ct * ( (dux(1)-duex(1))**2 )
         du2x = du2x + ct * ( (dux(2)-duex(2))**2 )
         du1y = du1y + ct * ( (duy(1)-duey(1))**2 )
         du2y = du2y + ct * ( (duy(2)-duey(2))**2 )
c
c    pressao descontinua
c
	 pe = 0.d00
	 pex = 0.d00
	 pey = 0.d00
	 do i=1,nenp
	  pe = pe + shgp(3,i,l)*dlp(1,i)
	  pex= pex+ shgp(1,i,l)*dlp(1,i)
	  pey= pey+ shgp(2,i,l)*dlp(1,i)
	 end do
	  dpe = dpe  + ct*((pe - ue(3))**2)
	  dpex= dpex + ct*((pex-duex(3))**2)
	  dpey= dpey + ct*((pey-duey(3))**2)
c
      if(iwrite.ne.0) then
c         write(43,2424) n,l,u(1),ue(1),u(2),ue(2),pe,ue(3)
2424   format(2i10,6e15.5)
222    format(8e20.9)
       end if

c
c     conservacao
c
      divue = divue + ct*(dux(1)-duex(1))**2 + ct*(duy(2)-duey(2))**2
c
4040  continue
c
c
	  edu  = edu  + (un1 + un2)
	  edux = edux + (du1x + du2x)
	  eduy = eduy + (du1y + du2y)
	  divu = divu + divue
	  edp  = edp  + dpe
	  edpx = edpx + dpex
	  edpy = edpy + dpey
c
c
c......boundary terms - multiplier
c
      xmlte = 0.d00
c
      do 4000 ns=1,nside
c
c    Material do lado
c
      ndgs = lado(ns,n) ! numero global do lado
c
c-----localiza os parametros do lado ns
c
      do nn=1,npars
         nld = (ns-1)*npars + nn
         do jj=1,ndof
	  dls(jj,nn) = dl(jj,nld)
         end do
      end do
c
c-----localiza os no's do lado ns
c
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,n)
      nl2=ien(ns2,n)
c      
      if(nl2.gt.nl1) then 
	  sign = 1.d00
        do  nn=1,nenlad
           idlsd(nn)=idside(ns,nn)
	  end do
	else
	  sign = -1.d00
	  idlsd(1) = idside(ns,2)
	  idlsd(2) = idside(ns,1)
	  id3 = nenlad-2
	  if(id3.gt.0) then
	    do il=3,nenlad
	      idlsd(il) = idside(ns,nenlad+3-il)
	    end do 
	  end if
      end if 
c
c
      do 2000 nn=1,nenlad
      nl=idlsd(nn)
      xls(1,nn)=xl(1,nl)
      xls(2,nn)=xl(2,nl)
c
 2000 continue
c
c
      call oneshgp(xls,detn,shlb,shln,shgn,
     &             nenlad,nnods,nints,nesd,ns,n,neg)      
c
      call oneshgp(xls,detpn,shlb,shlpn,shgpn,
     &             nenlad,npars,nints,nesd,ns,n,neg)      
c  
c
c.....compute boundary integral
c

      do 1000 ls=1,nints
c
      cwn = wn(ls)*detn(ls)
c
c    valores dos parametros do multiplicador
c
        dhs = 0.d00
          do i=1,npars
	      dhs = dhs + dls(1,i)*shgpn(2,i,ls)
          end do
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c
c    valores exatos dos multiplicadores
c

          call uexafx(x1,x2,ue,duex,duey,alpha,mm,index)
c
       dhse=ue(3)
c
c
      xmlte = xmlte + cwn*(dhse-dhs)**2
c
 1000 continue
 4000 continue
c
      xmlt = xmlt + xmlte
c
  50   continue
c
c




c...... error in all domain
c
          
       edu  = dlog10(dsqrt(edu))
       gradu= dlog10(dsqrt(edux+eduy))
       divu = dlog10(dsqrt(divu))
	   xmlt  = dlog10(dsqrt(xmlt))
c
	   edp = dlog10(dsqrt(edp))
	   gradp= dlog10(dsqrt(edpx+edpy))	   

c	 
       xel = numel
       xel = dlog10(xel)/2.d00
       xeq = neq
       xeq = dlog10(xeq)
c
      write(iplt,2101) xel,xeq,edu,gradu,edp,gradp,divu,xmlt
       return
c
2100  format(5x,i3,4x,e15.7,5x,e15.7)
2101  format(10(e13.4))
       end
c**** new **********************************************************************
c**** new **********************************************************************
      subroutine flninter(ien   ,x     ,xl   ,    
     &                 d     ,dl    ,mat   ,
     &                 c     ,ipar  ,dlf   ,   
     &                 dlp   ,dsfl  ,det   ,
     &                 shl   ,shg   ,wt    ,
     &                 detc  ,shlc  ,shgc  , 
     &                 ddis  ,detp  ,shlp  , 
     &                 shgp ,  
c
     &                 shln  ,shgn  ,
     &                 detn  ,shlb  ,shgb  ,
     &                 detpn ,shlpn ,shgpn ,
     &                 idside,xls   ,idlsd ,
     &                 grav  ,wn    ,
c                  
     &                 numel ,neesq ,nen   ,nsd   ,
     &                 nesd  ,nint  ,neg   ,nrowsh,
     &                 ned   ,nee   ,numnp ,ndof  ,
     &                 ncon  ,nencon,necon ,index ,
     &                 nints ,iwrite,iplt  ,ierb  ,
     &                 nenp   ,nside,nnods ,nenlad, 
     &                 npars  ,nmultp, nodsp)
c------------------------------------------------------------------
c     Some variables used in this subroutine
c            j: degree of freedom (1,...ncon)
c         u(j): finite element solution
c        du(j): derivative of finite element solution
c        ue(j): exact solution
c       due(j): derivative of exact solution
c       el2(j): error in L2
c      epri(j): error in the seminorm of H1 (L2 of derivatives)
c     el2el(j): error in L2 in the element domain
c    epriel(j): error in the seminorm of H1 in the element domain
c----------------------------------------------------------------
c
c
c    Program to calculate and print the L2 and H1 seminorm of the error
c    for each degree of freedom in the finite element solution .
c    The trapezoidal rule is used for integration in each element .
c    The number of integration points is given by nints .
c
c    This version is only applicable to 2D.
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      character *1 tab
      dimension ien(nen,*),x(nsd,*),xl(nesd,*),d(ndof,*),dl(ndof,*),
     &          mat(*),c(9,*),ipar(nodsp,*),dlf(ncon,*) ,
     &          dsfl(ncon,nencon,*),ddis(ned,nenp,*)
      dimension u(6),dux(6),duy(6),ue(6),duex(6),duey(6)
      dimension shl(3,nen,*),shg(3,nen,*),det(*),wt(*)
      dimension shlc(3,nencon,*),shgc(3,nencon,*),detc(*)
      dimension el2(6),epri(6),eprix(6),epriy(6),el2el(6),
     &          eprxel(6),epryel(6)
	  dimension detp(*),shlp(3,nenp,*),shgp(3,nenp,*),dlp(ned,*)
c
      dimension detn(*),shlb(3,nenlad,*),shgb(3,nenlad,*),
     &          detpn(*),shlpn(3,npars,*),shgpn(3,npars,*),
     &          idside(nside,*),xls(nesd,*),idlsd(*)
	  dimension dls(ndof,16),grav(*),wn(*)
      dimension shln(3,nnods,*),shgn(3,nnods,*)
      dimension sxlhp(64,64),sxlhv(64,64),sxlhs(8),xlpn(2,64),xlvn(2,64)
      dimension shlsd(8,8),xlps(2,8)
c
      common /colhtc/ neq
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
      pi=4.*datan(1.d00)
      dpi=2.d00*pi
c
      gf1=grav(1)
      gf2=grav(2)
      gf3=grav(3) 
      gf4=grav(4)
      gf5=grav(5) 
      gf6=grav(6)      
        tab=char(9)
         xint=0.d00 
         yint=0.d00
c
      call clear ( el2, ncon )
      call clear (eprix, ncon )
      call clear (epriy, ncon )
c
      divu=0.d00
	  edp = 0.d00
	  edpx= 0.d00
	  edpy= 0.d00
	  xmlt= 0.d00
c
c
       call shapesd(shlsd,nenlad,npars)
      if (nen.eq.3) then
       call shapent(sxlhp,nen,nenp)
       call shapent(sxlhv,nen,nencon)
      else if (nen.eq.4) then
        call shapen(sxlhp,nen,nenp)
        call shapen(sxlhv,nen,nencon)
      else
        stop
      endif
        

c
c.....  loop on elements
c      
      do 50 n=1,numel
c
        mm = mat(n)
        alpha=c(9,mm)
	    xka1=c(1,mm)*alpha  !permeabilidade 1
	    xka2=c(2,mm)*alpha  !permeabilidade 2
	    eps=c(3,mm)  !
c
        xkl = c(1,mm)
c
       call local(ien(1,n),x,xl,nen,nsd,nesd)
       call local(ipar(1,n),d,dl,nodsp,ndof,ndof)
c
      do i=1,nenp
      xlpn(1,i) = 0.d00
      xlpn(2,i) = 0.d00
      do j=1,nen
        xlpn(1,i) = xlpn(1,i) + sxlhp(j,i)*xl(1,j)
        xlpn(2,i) = xlpn(2,i) + sxlhp(j,i)*xl(2,j)
       end do
       end do
c       
       call elemdlp(xlpn,dlp,alpha,ned,nenp,mm,index)
c
      do i=1,nencon
      xlvn(1,i) = 0.d00
      xlvn(2,i) = 0.d00
      do j=1,nen
        xlvn(1,i) = xlvn(1,i) + sxlhv(j,i)*xl(1,j)
        xlvn(2,i) = xlvn(2,i) + sxlhv(j,i)*xl(2,j)
       end do
       end do
c       
       call elemdlf(xlvn,dlf,ncon,nencon,alpha,mm,index)
c
      call clear ( el2el, ncon )
      call clear (eprxel, ncon )
      call clear (epryel, ncon )
      divue=0.d00
	  pmede=0.d00
	  dpe = 0.d00
	  dpex= 0.d00
	  dpey= 0.d00
c
c.....  loop on integration points
c
c
c.....Triangles or quadrilaterals
c
c
      call shgqs(xl,detc,shlc,shgc,nint,n,neg,.true.,nencon,shl,nen)
c
      call shgqs(xl,detp,shlp,shgp,nint,n,neg,.true.,nenp,shl,nen)
c
c
        do 4040 l=1,nint
        ct=detc(l)*wt(l)
        call clear ( u, ncon )
        call clear (dux, ncon )
        call clear (duy, ncon )
        xint=0.d0
        yint=0.d0
         do 3030 i=1,nen
            xint=xint+shl(3,i,l)*xl(1,i)
            yint=yint+shl(3,i,l)*xl(2,i)
3030     continue
c
          do 3330 i=1,nencon 
            do 2020 j=1,ncon
                u(j)   = u(j) + shgc(3,i,l)*dlf(j,i)
	          dux(j) = dux(j) + shgc(1,i,l)*dlf(j,i)
	          duy(j) = duy(j) + shgc(2,i,l)*dlf(j,i)
2020        continue
3330      continue
c

      call uexafx(xint,yint,ue,duex,duey,alpha,mm,index)
c
c    pressao descontinua
c
	pe = 0.d00
	pex = 0.d00
	pey = 0.d00
	do i=1,nenp
	  pe = pe + shgp(3,i,l)*dlp(1,i)
	  pex= pex+ shgp(1,i,l)*dlp(1,i)
	  pey= pey+ shgp(2,i,l)*dlp(1,i)
	end do
	  dpe = dpe + ct*(pe - ue(3))**2
	  dpex= dpex + ct*(pex-duex(3))**2
	  dpey= dpey + ct*(pey-duey(3))**2
c
      if(iwrite.ne.0) then
c         write(43,2424) n,l,u(1),ue(1),u(2),ue(2),pe,ue(3)
2424   format(2i10,6e15.5)
222    format(8e20.9)
       end if
c
       do 3535 j=1,ncon
        un = ct * ( (u(j)-ue(j))**2 )
        upnx= ct * ( (dux(j)-duex(j))**2 )
        upny= ct * ( (duy(j)-duey(j))**2 )
        el2el(j) = el2el(j) + un
        eprxel(j) = eprxel(j) + upnx
        epryel(j) = epryel(j) + upny
3535   continue
c
c     conservacao
c
      divue = divue + ct*(dux(1)-duex(1))**2 + ct*(duy(2)-duey(2))**2
c
4040  continue
c
       do 45 j=1,ncon
        el2(j) = el2(j) + el2el(j)
        eprix(j) = eprix(j) + eprxel(j)
        epriy(j) = epriy(j) + epryel(j)
  45   continue
c
	   divu = divu + divue
	   edp = edp + dpe
 	   edpx= edpx+ dpex
	   edpy= edpy+ dpey
c
c......boundary terms - multiplier
c
      xmlte = 0.d00
c
      do 4000 ns=1,nside
c
c-----localiza os parametros do lado n
c
      do nn=1,npars
         nld = (ns-1)*npars + nn
         do jj=1,ndof
	  dls(jj,nn) = dl(jj,nld)
         end do
      end do
c
c-----localiza os no's do lado ns
c
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,n)
      nl2=ien(ns2,n)
c      
      if(nl2.gt.nl1) then 
	  sign = 1.d00
        do  nn=1,nenlad
           idlsd(nn)=idside(ns,nn)
	  end do
	else
	  sign = -1.d00
	  idlsd(1) = idside(ns,2)
	  idlsd(2) = idside(ns,1)
	  id3 = nenlad-2
	  if(id3.gt.0) then
	    do il=3,nenlad
	      idlsd(il) = idside(ns,nenlad+3-il)
	    end do 
	  end if
      end if 
c
c
      do 2000 nn=1,nenlad
      nl=idlsd(nn)
      xls(1,nn)=xl(1,nl)
      xls(2,nn)=xl(2,nl)
c
 2000 continue
c
      do nn=1,npars
      xlps(1,nn) = 0.d00
      xlps(2,nn) = 0.d00
      do jj=1,nenlad
        xlps(1,nn) = xlps(1,nn) + shlsd(jj,nn)*xls(1,jj)
        xlps(2,nn) = xlps(2,nn) + shlsd(jj,nn)*xls(2,jj)
      end do
      end do
c
       call elemdls(xlps,dls,ndof,npars,alpha,mm,index)
c
      call oneshgp(xls,detn,shlb,shln,shgn,
     &             nenlad,nnods,nints,nesd,ns,n,neg)      
c
      call oneshgp(xls,detpn,shlb,shlpn,shgpn,
     &             nenlad,npars,nints,nesd,ns,n,neg)      
c  
c
c.....compute boundary integral
c

      do 1000 ls=1,nints
c
      cwn = wn(ls)*detn(ls)
c
c    valores dos parametros do multiplicador
c
        dhs = 0.d00
          do i=1,npars
	      dhs = dhs + dls(1,i)*shgpn(2,i,ls)
          end do
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c
c    valores exatos dos multiplicadores
c

         call uexafx(x1,x2,ue,duex,duey,alpha,mm,index)
c
c
c
       dhse=ue(3)
c
      xmlte = xmlte + cwn*(dhse-dhs)**2
c
 1000 continue
 4000 continue
c
      xmlt = xmlt + xmlte
c
  50   continue
c
c
c...... error in all domain
c
       fxl2 = dsqrt(el2(1) + el2(2))
       fxh1 = dsqrt(eprix(1)+eprix(2) + epriy(1)+epriy(2))
      
	   edp = dlog10(dsqrt(edp))
	   gradp= dlog10(dsqrt(edpx+edpy))

       fxl2 = dlog10(fxl2)
       fxh1 = dlog10(fxh1)
	   xmlt = dlog10(dsqrt(xmlt))
       xel=numel
       if(nen.eq.3.or.nen.eq.6.or.nen.eq.10) xel=numel/2.
       xel = dlog10(xel)/2.d00
       xeq = neq
       xeq = dlog10(xeq)
c
       fx1l2 = dsqrt(el2(1))
       fx2l2 = dsqrt(el2(2))
       fx1l2 = dlog10(fx1l2)
       fx2l2 = dlog10(fx2l2)
    	divu = dlog10(dsqrt(divu))

      write(iplt,2101) xel,xeq,fxl2,fxh1,edp,gradp,divu,xmlt
       return
c
2100  format(5x,i3,4x,e15.7,5x,e15.7)
2101  format(12(e13.4))
       end
c**** new **********************************************************************
c**** new **********************************************************************
c-----------------------------------------------------
      subroutine elemdlf(xl,dlf,ncon,nencon,alpha,m,index)
c-----------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension xl(2,*),dlf(ncon,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      pi=4.d00*datan(1.d00)
      pi2=pi*pi
c
      do 1001 n=1,nencon 
       x = xl(1,n)
       y = xl(2,n)
c
       go to (100,200,300) index
c
c
c        PROBLEM 1 - CILAMCE 2014 (CRUMPTON)
c
 100   continue
c
      px=pi*x
      py=pi*y
      sx=dsin(px)
      sy=dsin(py)
      cx=dcos(px)
      cy=dcos(py)
      ex=dexp(x)
c
      if(m.eq.1) then
c
c       Solution Homogeneous Darcy
c
c
      pia2 = 2.d00*pi
c     u1 e u2      
      dlf(1,n) = -pia2*cx*sy
      dlf(2,n) = -pia2*sx*cy
c
      else
c
c       Solution anisotropic Darcy
c
c     u1 e u2      
      pia = pi*alpha
      dlf(1,n) = -pia*(2.d00*cx*sy + sx*cy)
      dlf(2,n) = -pia*(cx*sy + 2.d00*sx*cy)
c
      end if
c
      go to 1001
c
c        PROBLEM 2 - CILAMCE 2014 (CRUMPTON)
c
 200   continue
c
      sx=dsin(x)
      sy=dsin(y)
      cx=dcos(x)
      cy=dcos(y)
      ex=dexp(x)
c
      if(m.eq.1) then
c
      dlf(1,n) = -alpha*(2.d00*sy+cy)
      dlf(2,n) = -alpha*x*(2.d00*cy-sy)-cy
c
      else
c
c       Solution anisotropic Darcy
c
c     u1 e u2      
      dlf(1,n) = -alpha*(2.d00*ex*sy + ex*cy)
      dlf(2,n) = -alpha*(ex*sy + 2.d00*ex*cy)
c
      end if
      go to 1001
c
 300   continue
c
      px=pi*x
      py=pi*y
      sx=dsin(px)
      sy=dsin(py)
      cx=dcos(px)
      cy=dcos(py)
      pi2=pi*pi
      co=1.d00/(pi2*2.d00)
c      
      dlf(1,n) = co*pi*sx*cy
      dlf(2,n) = co*pi*cx*sy
c
1001  continue
      return
      end
c**** new ********************************************************************** 

c-----------------------------------------------------
      subroutine elemdlp(xl,dlp,alpha,ned,nenp,m,index)
c-----------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension xl(2,*),dlp(ned,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      pi=4.d00*datan(1.d00)
      pi2=pi*pi
c
      do 1001 n=1,nenp 
       x = xl(1,n)
       y = xl(2,n)
c
       go to (100,200,300) index
c
c
c        PROBLEM 1 - CILAMCE 2014 (CRUMPTON)
c
 100   continue
c
      px=pi*x
      py=pi*y
      sx=dsin(px)
      sy=dsin(py)
      cx=dcos(px)
      cy=dcos(py)
      ex=dexp(x)
c
      if(m.eq.1) then
c
c       Solution Homogeneous Darcy
c
c     p
c
      dlp(1,n) = 2.d00*sx*sy
c
c
      else
c
c       Solution anisotropic Darcy
c
c     p
c
      dlp(1,n) = sx*sy

c
      end if
c
      go to 1001
c
c        PROBLEM 2 - CILAMCE 2014 (CRUMPTON)
c
 200   continue
c
      sx=dsin(x)
      sy=dsin(y)
      cx=dcos(x)
      cy=dcos(y)
      ex=dexp(x)

      if(m.eq.1) then
c
c       Solution Homogeneous Darcy
c
c     p
c
      dlp(1,n) = alpha*x*(2.d00*sy+cy)+sy
c
c
      else
c
c       Solution anisotropic Darcy
c
c     p
c
      dlp(1,n) = ex*sy

c
      end if
      go to 1001
 300   continue
c
      px=pi*x
      py=pi*y
      sx=dsin(px)
      sy=dsin(py)
      cx=dcos(px)
      cy=dcos(py)
      pi2=pi*pi
      co=1.d00/(pi2*2.d00)
c      
      dlp(1,n) = co*cx*cy
c
1001  continue
      return
      end
c**** new **********************************************************************
c-----------------------------------------------------
      subroutine elemdls(xl,dls,ndof,npars,alpha,m,index)
c-----------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension xl(2,*),dls(ndof,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      pi=4.d00*datan(1.d00)
      pi2=pi*pi
c
      do 1001 n=1,npars
       x = xl(1,n)
       y = xl(2,n)
c
       go to (100,200,300) index
c
c
c        PROBLEM 1 - CILAMCE 2014 (CRUMPTON)
c
 100   continue
c
      px=pi*x
      py=pi*y
      sx=dsin(px)
      sy=dsin(py)
      cx=dcos(px)
      cy=dcos(py)
      ex=dexp(x)

      if(m.eq.1) then
c
c       Solution Homogeneous Darcy
c
c
      pia2 = 2.d00*pi
c     u1 e u2      
      dls(1,n) = 2.d00*sx*sy
c
      else
c
c       Solution anisotropic Darcy
c
c     u1 e u2      
      pia = pi*alpha
      dls(1,n) = sx*sy
c
      end if
c
      go to 1001
c
c        PROBLEM 2 - CILAMCE 2014 (CRUMPTON)
c
 200   continue
c
      sx=dsin(x)
      sy=dsin(y)
      cx=dcos(x)
      cy=dcos(y)
      ex=dexp(x)
c
      if(m.eq.1) then
c
      dls(1,n) = alpha*x*(2.d00*sy+cy)+sy
c
      else
c
c       Solution anisotropic Darcy
c
c     u1 e u2      
      dls(1,n) = ex*sy
c
      end if
      go to 1001
c
 300   continue
c
      px=pi*x
      py=pi*y
      sx=dsin(px)
      sy=dsin(py)
      cx=dcos(px)
      cy=dcos(py)
      pi2=pi*pi
      co=1.d00/(pi2*2.d00)
c      
      dls(1,n) = co*cx*cy
c
1001  continue
      return
      end
c**** new **********************************************************************
c-----------------------------------------------------
      subroutine uexafx(x,y,ue,duex,duey,alpha,m,index)
c-----------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension ue(*),duex(*),duey(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
c
       go to (100,200,300) index
c
c
c        PROBLEM 1 - CILAMCE 2014 (CRUMPTON)
c
 100   continue
c
      pi=4.d00*datan(1.d00)
      pi2=pi*pi
      px=pi*x
      py=pi*y
      sx=dsin(px)
      sy=dsin(py)
      cx=dcos(px)
      cy=dcos(py)
      ex=dexp(x)

      if(m.eq.1) then
c
c       Solution Homogeneous Darcy
c
c
      ue(1) = -2.d00*pi*cx*sy
      ue(2) = -2.d00*pi*sx*cy
      ue(3) =  2.d00*sx*sy
c      
      duex(1) =  2.d00*pi2*sx*sy
      duex(2) = -2.d00*pi2*cx*cy
      duex(3) =  2.d00*pi*cx*sy
c
      duey(1) = -2.d00*pi2*cx*cy
      duey(2) =  2.d00*pi2*sx*sy
      duey(3) =  2.d00*pi*sx*cy
c  
c
      else
c
c       Solution anisotropic Darcy
c
      pia = pi*alpha
      ue(1) = -pia*(2.d00*cx*sy + sx*cy)
      ue(2) = -pia*(cx*sy + 2.d00*sx*cy)
      ue(3) =  sx*sy
c      
      duex(1) = -pia*pi*(-2.d00*sx*sy + cx*cy)
      duex(2) = -pia*pi*(-sx*sy + 2.d00*cx*cy)
      duex(3) =  pi*cx*sy
c
      duey(1) = -pia*pi*(2.d00*cx*cy - sx*sy)
      duey(2) = -pia*pi*(cx*cy - 2.d00*sx*sy)
      duey(3) =  pi*sx*cy  
      end if
c
      go to 1001
c
c        PROBLEM 2 - CILAMCE 2014 (CRUMPTON)
c
 200   continue
c
      sx=dsin(x)
      sy=dsin(y)
      cx=dcos(x)
      cy=dcos(y)
      ex=dexp(x)

      if(m.eq.1) then
c
      ue(1) = -alpha*(2.d00*sy+cy)
      ue(2) = -alpha*x*(2.d00*cy-sy)-cy
      ue(3) =  alpha*x*(2.d00*sy+cy)+sy
c
      duex(1) =  0.d00
      duex(2) = -alpha*(2.d00*cy-sy)
      duex(3) =  alpha*(2.d00*sy+cy)
c
      duey(1) = -alpha*(2.d00*cy-sy)
      duey(2) =  alpha*x*(2.d00*sy+cy)+sy
      duey(3) =  alpha*x*(2.d00*cy-sy)+cy
c  
c
      else
c
      ue(1) = -alpha*(2.d00*ex*sy + ex*cy)
      ue(2) = -alpha*(ex*sy + 2.d00*ex*cy)
      ue(3) =  ex*sy
c      
      duex(1) = -alpha*(2.d00*ex*sy + ex*cy)
      duex(2) = -alpha*(ex*sy + 2.d00*ex*cy)
      duex(3) =  ex*sy
c
      duey(1) = -alpha*(2.d00*ex*cy - ex*sy)
      duey(2) = -alpha*(ex*cy - 2.d00*ex*sy)
      duey(3) =  ex*cy
      end if
      go to 1001
 300   continue
c
      pi=4.d00*datan(1.d00)
      pi2=pi*pi
      px=pi*x
      py=pi*y
      sx=dsin(px)
      sy=dsin(py)
      cx=dcos(px)
      cy=dcos(py)
      co=1.d00/(pi2*2.d00)
c      
      ue(1) = co*pi*sx*cy
      ue(2) = co*pi*cx*sy
      ue(3) = co*cx*cy
c      
      duex(1)= co*pi2*cx*cy
      duex(2)=-co*pi2*sx*sy
      duex(3)=-co*pi*sx*cy
c      
      duey(1)=-co*pi2*sx*sy
      duey(2)= co*pi2*cx*cy
      duey(3)=-co*pi*cx*sy
c
1001  continue
      return
      end
	  
c***********************************************************************
c          impressões
c***********************************************************************
c**** new **********************************************************************
      subroutine printd(name,dva,ndof,numnp,icode)
c
c
c.... program to print kinematic data
c
      implicit real*8(a-h,o-z)
c
c.... remove above card for single precision operation
c
      logical lzero
	character *4 name(11)
      dimension dva(ndof,*)
c
      nn = 0
c

      do 100 n=1,numnp
      call ztest(dva(1,n),ndof,lzero)
      if (.not.lzero) then
         nn = nn + 1
         if (mod(nn,50).eq.1) 
     &      write(icode,1000) name,(i,i=1,ndof)
         write(icode,2000) n,(dva(i,n),i=1,ndof)
      endif
  100 continue
c
      return
c
 1000 format(///,11a4//6xx,'node',6(11x,'dof',i1)/)
 2000 format(1x,i10,2x,6(1pe13.6,2x))
      end
c**** new **********************************************************************
      subroutine printdX(name,dva,ndof,numnp,icode)
c     Igual a printd (Boness) criei para fazer testes
c
c.... program to print kinematic data
c
      implicit real*8(a-h,o-z)
c
c.... remove above card for single precision operation
c
      logical lzero
	character *4 name(11)
      dimension dva(ndof,*)
c
      nn = 0
c

      do 100 n=1,numnp
c     call ztest(dva(1,n),ndof,lzero)
c     if (.not.lzero) then
         nn = nn + 1
!          if (mod(nn,50).eq.1) 
!      &      write(icode,1000) name,(i,i=1,ndof)
         write(icode,2000) n,(dva(i,n),i=1,ndof)
c     endif
  100 continue
c
      return
c
 1000 format(///,11a4//6xx,'node',6(11x,'dof',i1)/)
 2000 format(1x,i10,2x,6(1pe13.6,2x))
      end
      
c***********************************************************************
      subroutine printf(f,ndof,numnp,nlv)
c
c.... program to print prescribed force and boundary condition data
c
      implicit real*8(a-h,o-z)
c
c.... remove above card for single precision operation
c
      logical lzero
      dimension f(ndof,numnp,*)
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
      nn = 0
c
      do 100 n=1,numnp
      call ztest(f(1,n,nlv),ndof,lzero)
      if (.not.lzero) then
         nn = nn + 1
         if (mod(nn,50).eq.1) 
     &      write(iecho,1000) nlv,(i,i=1,ndof)
         write(iecho,2000) n,(f(i,n,nlv),i=1,ndof)
      endif
  100 continue
c
      return
c
 1000 format(///,
     &' p r e s c r i b e d   f o r c e s   a n d   k i n e m a t i c ',
     &'  b o u n d a r y   c o n d i t i o n s'//5x,
     &' load vector number = ',i10///5x,
     &' node no.',6(13x,'dof',i1,:)/)
 2000 format(6x,i10,10x,6(1pe15.8,2x))
      end
c**** new ********************************************************************** 
      subroutine prntel(mat,ien,nen,numel)
c
c.... program to print data for element with "nen" nodes
c
c        note: presently the label formats are limited to
c              elements with one to nine nodes
c
      dimension mat(*),ien(nen,*)
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
      do 100 n=1,numel
      if (mod(n,50).eq.1) write(iecho,1000) (i,i=1,nen)
      write(iecho,2000) n,mat(n),(ien(i,n),i=1,nen)
  100 continue
c
      return
c
 1000 format(///,
     &' d a t a   f o r   e l e m e n t   n o d e s ',//1x,
     &' element      material',16('      node',i2))
 2000 format(1x,i10,20(2x,i10))
      end
c
c**** new **********************************************************************
      subroutine matlado(lado,mat,matside,mcont,nside,numel,nedge)
c
c.... program to print data for element with "nen" nodes
c
c        note: presently the label formats are limited to
c              elements with one to nine nodes
c
      dimension lado(nside,*),mat(*),matside(*),mataux(nedge),mcont(*)
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
      do 100 n=1,numel
       do ns=1,nside
         na = lado(ns,n)
         mataux(na) = mat(n)
         mcont(na)=0
       end do
       do ns=1,nside
         na = lado(ns,n)
       if(mataux(na).gt.matside(na)) then
         matside(na) = mataux(na)
       end if
       end do
  100 continue
       do n=1,numel
       do ns=1,nside
         na = lado(ns,n)
       if(mat(n).ne.matside(na))then
       mcont(na)=1
       end if
       end do
       end do
      write(iecho,1000) nedge
      do i=1,nedge
       write(iecho,2000) i, matside(i)
      end do
c
      return
c
 1000 format(///, i10,
     &' m a t e r i a l  g r u p   o f  s i d e s ',//1x,
     &' lado      material grup')
 2000 format(1x,2i10)
      end
c**** new **********************************************************************
      subroutine prntels(mat,ien,nen,numel)
c
c.... program to print data for element with "nen" nodes
c
c        note: presently the label formats are limited to
c              elements with one to nine nodes
c
      dimension mat(*),ien(nen,*)
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
      do 100 n=1,numel
      if (mod(n,50).eq.1) write(iecho,1000) (i,i=1,nen)
      write(iecho,2000) n,mat(n),(ien(i,n),i=1,nen)
  100 continue
c
      return
c
 1000 format(///,
     &' d a t a   f o r    e l e m e n t   s i d e s ',//1x,
     &' element      material',16('      node',i2))
 2000 format(1x,i10,20(2x,i10))
      end
c**** new **********************************************************************
      subroutine prntelp(mat,ien,nen,numel)
c
c.... program to print data for element with "nen" nodes
c
c        note: presently the label formats are limited to
c              elements with one to nine nodes
c
      dimension mat(*),ien(nen,*)
      common /iounit/ iin,iout,iecho,ioupp,iout1,itest2,ierrb0,ierrb,
     & ierrbi
c
      do 100 n=1,numel
      if (mod(n,50).eq.1) then
	  if(nen.le.16) write(iecho,1000) (i,i=1,nen)
	  if(nen.eq.20) write(iecho,1020) (i,i=1,nen)
	  if(nen.eq.24) write(iecho,1024) (i,i=1,nen)
	  if(nen.eq.28) write(iecho,1028) (i,i=1,nen)
	  if(nen.eq.32) write(iecho,1032) (i,i=1,nen)
	end if
      write(iecho,2000) n,mat(n),(ien(i,n),i=1,nen)
  100 continue
c
      return
c
 1000 format(///,
     &' d a t a   f o r   e l e m e n t   p a r a m e t e r s ',//1x,
     &' element      material',16('      node',i2))
 1020 format(///,
     &' d a t a   f o r   e l e m e n t   p a r a m e t e r s ',//1x,
     &' element      material',20('      node',i2))
 1024 format(///,
     &' d a t a   f o r   e l e m e n t   p a r a m e t e r s ',//1x,
     &' element      material',24('      node',i2))
 1028 format(///,
     &' d a t a   f o r   e l e m e n t   p a r a m e t e r s ',//1x,
     &' element      material',28('      node',i2))
 1032 format(///,
     &' d a t a   f o r   e l e m e n t   p a r a m e t e r s ',//1x,
     &' element      material',32('      node',i2))
2000  format(1x,i10,64(2x,i10))
      end
c**** new **********************************************************************
      subroutine printdlocal(name,dva,ndof,numnp,numel,icode)
c
c
c.... program to print kinematic data
c
      implicit real*8(a-h,o-z)
c
c.... remove above card for single precision operation
c
      logical lzero
	character *4 name(11)
      dimension dva(ndof,numnp,*)
c     
c      write(icode,1001) nel
c      write(icode,1000) name,(i,i=1,ndof)
        
      do nel=1,numel
         nn = 0
         do n=1,numnp
            nn = nn + 1
             if (mod(nn,numel).eq.1)  then
              write(icode,1001) nel
              write(icode,1000) name,(i,i=1,ndof)
             endif 
            write(icode,2000) n,(dva(i,n,nel),i=1,ndof)            
         end do
      end do
c     
      return
c
 1000 format(11a4//6xx,'node',6(11x,'dof',i1)/)
 1001 format(//,6xx,'element',i5)
 2000 format(1x,i10,2x,6(1pe13.6,2x))
  
      end
c**** new **********************************************************************
c**** PLOT**********************************************************************
c**** new **********************************************************************
      subroutine aproxprn(dsfl, ddis, x, ien, mat,
     1                    ncon, nencon, numel, nenp,
     2                    nsd, nen, numnp, neg, numat,
     3                    nesd, iecho,c)

      implicit real*8(a-h,o-z)

      !!! ATEN‚ÌO: imprime apenas elementos lineares e bilineares!!!
      ! dsfl(ncon, nencon, numel): local velocity field
      ! ddis(nenp, numel): local pressure field

	  dimension ien(nen, numel),c(9,*)
      dimension dsfl(ncon, nencon, numel)
      dimension ddis(nenp, numel)
      dimension x(nsd, numnp)
      dimension xl(nesd, nen),ue(6),duex(6),duey(6)
      do i = 1, numel
         do j = 1, nen
            k = ien(j, i)
            do m=1,nsd
               xl(m, j) = x(m, k)
            end do
         end do
		 
         do j = 1, nen
            k = ien(j, i)
c
         x1=x(1, k)
         x2=x(2, k)
c         call uexafx(x1,x2,ue,duex,duey,eps,stk,drc)
c
            write(2, "(8e15.7)") x(1, k), x(2, k),
     1                           ue(1),ue(2),ue(3)

c
            write(1, "(8e15.7)") x(1, k), x(2, k),
     1                             (dsfl(ndof, j, i), ndof=1,2),
     2                             ddis(j, i)
         end do
c         ! reimprime o primeiro no para fechar um plano no plot
         j = 1
         k = ien(j, i)
c
         x1=x(1, k)
         x2=x(2, k)
c         call uexafx(x1,x2,ue,duex,duey,eps,stk,drc)
c
            write(2, "(8e15.7)") x(1, k), x(2, k),
     1                           ue(1),ue(2),ue(3)


         write(1, "(8e15.7)") x(1, k), x(2, k),
     1                          (dsfl(ndof, j, i), ndof=1,2),
     2                          ddis(j, i)
         write(1, "(/)")
         write(2, "(/)")
      end do

      write(1, "(/)")
      write(1, "(/)")
      write(2, "(/)")
      write(2, "(/)")

      


      return
      end
