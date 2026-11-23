      subroutine rosame_filtre_demerliac(nb_val_max,nb_val,rtab_val_in,
     -                                   rflag,rtab_val_out)
!     -----------------------------------------------------------------

! filtrage d'une serie temporelle par un filtre de demerliac a 3 jours pour
! enlever tout ce qui est inferieur a 3 jours
! ne fonctionne que pour des donnees horaires
! creation Philippe Techine (version 0.0 - juillet 1999)
! a partir de demerliac cree par Jean-Marc Molines (1998)
!
! Routine for 3 days Demerliac filter for hourly time series.
! Demerliac filter removes all scale lower than 3 days on time series .
! Only works on hourly time series .
! written by Philippe Techine (version 0.0 - juillet 1999)
! and Jean-Marc Molines (1998)

      implicit none

! variables entree
! input variables 
      integer nb_val_max,nb_val
      real rtab_val_in(nb_val_max)
      real rflag

! variables sortie
! output variables
      real rtab_val_out(nb_val_max)

! variables locales
! local variables 
      integer i,k
      integer idebut,ifin

! 24 heures * 3 jours -1 = 71 points
! 24 hours * 3 days -1 = 71 points
      integer nbre
      data nbre/71/
      integer itab(71),idata(36)
      data idata/1,3,8,15,21,32,45,55,72,91,105,128,153,171,
     -           200,231,253,288,325,351,392,435,465,512,558,
     -           586,624,658,678,704,726,738,752,762,766,768/
      real somme

! initialisations
      k=0
      somme=0.
      idebut=(nbre-1)/2+1
      ifin=nb_val-idebut+1

! remplissage du tableau et calcul de la somme
! fill the array and compute sum 
      do i=1,(nbre+1)/2
         k=k+1
         itab(k)=idata(i)
      enddo
      do i=(nbre-1)/2,1,-1
         k=k+1
         itab(k)=idata(i)
      enddo
      do i=1,nbre
         somme=somme+real(itab(i))
      enddo

! filtrage des valeurs
! data filtering
      do i=idebut,ifin
         rtab_val_out(i)=0.
         do k=1,nbre
            rtab_val_out(i)=rtab_val_out(i)
     -                     +real(itab(k))*rtab_val_in(i-idebut+k)/somme
         enddo
      enddo

! les 35 premieres et 35 dernieres valeurs sont flaguees pour eviter les  effets de bord
! 35 first and last values of time series are flagged
      do i=1,idebut-1
         rtab_val_out(i)=rflag
         rtab_val_out(nb_val-i+1)=rflag
      enddo
      end
