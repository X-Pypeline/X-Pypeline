# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017)
#
# This file is part of XPypeline.
#
# GWpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GWpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GWpy.  If not, see <http://www.gnu.org/licenses/>.

"""This module converts ra to dec and saves to a file if supplied
"""

import lal
import numpy as np

def radectoearth(ra, dec, gps, filename=''):
    gmst = lal.GreenwichMeanSiderealTime(gps)
    gmst_deg = gmst / 86400 * 360;
    # ---- Compute Earth-based coordinates of targeted sky location, in degrees.
    phi_deg = ra - gmst_deg;
    theta_deg  = 90 - dec;
    # ---- Convert to radians.
    phi = phi_deg * 2 * np.pi / 360;
    phi = np.mod(phi+np.pi,2*np.pi)-np.pi;
    theta = theta_deg * 2 * np.pi / 360;
    if filename:
        with open(filename, 'w') as f:
            f.write('{0} {1}'.format(phi, theta))
            f.close()
        return
    else:
        return phi, theta


def xmakeskygrid(ra_ctr_deg, dec_ctr_deg, gps, 
                 sigma_deg, nSigma, sites, delay_tol, outputfile=None, gridtype='circular', verbose=False):
    """
    # XMAKESKYGRID - tile the sky to cover circular error boxes
    #
    # XMAKESKYGRID - Given a set of estimated sky locations and
    # common error in sky location measurement, generate a list of sky
    # locations which we should search over to keep time-delay errors below a
    # given threshold for the specified network.
    #
    # usage:
    #
    # [ra_search_deg,dec_search_deg,probabilities,gridarea] = ...
    #     xmakeskygrid(ra_ctr_deg,dec_ctr_deg,gps,sigma_deg, ...
    #     nSigma,sites,delay_tol,outputfile,gridtype,verbose)
    #
    #  ra_ctr_deg        Tilde-delimited string. Right ascension of centre of 
    #                    error circles [deg], e.g., '178.1~172.3'.
    #  dec_ctr_deg       Tilde-delimited string. Declination of centre of
    #                    error circles [deg], e.g., '-12.0~-22.2'.
    #  gps               String. GPS time of trigger.
    #  sigma_deg         Tilde-delimited string. 1-sigma size of error circles
    #                    [deg]. This is
    #                    multiplied by nSigma (see below) to define a radius
    #                    around each central point to be covered by the grid.
    #                    A single value may be supplied, in which case this
    #                    common value is used for all error circles.
    #  nSigma            String. Determines region about each central position
    #                    to be covered by grid.
    #  sites             String. Tilde-separated list of detector sites,
    #                    e.g., 'H~L~V'.
    #  delay_tol         String. Maximum error in time delay [sec] allowed
    #                    between any sky position in an error box and the
    #                    closest grid point. 
    #  outputfile        String, optional.  Name of output file to which to
    #                    write grid. The file will contain four columns:
    #                    1. theta   Polar angle [rad] in the range [0, pi] with
    #                       0 at the North Pole/ z axis.
    #                    2. phi     Azimuthal angle [rad]in the range [-pi, pi) 
    #                       with 0 at Greenwich / x axis.
    #                    3. pOmega  A priori probability associated with each
    #                       grid point.
    #                    4. dOmega  Solid angle associated with each grid point.
    #                    This follows the format of the skyPositions array as
    #                    described in sinusoidalMap.m. Specify 'None' for no
    #                    file to be produced.
    #  gridtype          String, optional.  Type of grid.  Recognized values:
    #                      'circular'       Grid made of concentric rings centered
    #                        on each input sky position, with redundant points
    #                        discarded.  Made using SINUSOIDALMAP.
    #                      'healpix'    Single uniform all-sky grid, with
    #                        points outside all error circles discarded.  Made
    #                        using HEALPIX.
    #  verbose           String, optional.  Use '1' for verbose output.
    #                    Default '0'.
    #
    #  ra_search_deg     Vector. Right ascensions [deg] of search grid points.   
    #  dec_search_deg    Vector. Declinations [deg] of search grid points.
    #  probabilities     Vector. Fisher probability of each sky position. For
    #                    1-sigma errors above 360 degrees uniform probability
    #                    is used.
    #  gridarea          Vector. Area [steradians] of each grid point.
    #
    #  The probability for each grid point is the sum of the Fisher
    #  probabilities for that grid point with each of the input central sky
    #  positions. 
    #
    # $Id: xmakeskygrid.m 4364 2014-07-19 07:20:42Z iain.dorrington@LIGO.ORG $
    """
    """
    ###########################################################################
    #                              Preliminaries.
    ###########################################################################
    # ---- Convert input variables to appropriate type.
    sites = sites.split('~')

    # ---- Length checks.
    if len(ra_ctr_deg) != len(dec_ctr_deg):
        raise ValueError('Right ascension and declination must have same length.')

    if len(sigma_deg) == 1:
        sigma_deg = sigma_deg * np.ones(len(ra_ctr_deg))

    elif len(sigma_deg) ~= len(ra_ctr_deg):
        raise ValueError('Sigma must be scalar or same length as right ascension, declination.')

    # ---- Speed of light (in m/s).
    speedOfLight = 299792458;
    # ---- Calc radius of sky area we will search for GRBs.
    skyPosError_deg = sigma_deg * nSigma;
    # ---- Radius of sky area in radians.
    skyPosError_rad = skyPosError_deg * pi/180;
    ###########################################################################
    #          Get cartesian Earth-based coordinates (m) of detector
    #                       vertex for each site.
    ###########################################################################
    nSites = len(sites)
    for iSite = range(nSites):
        det = LoadDetectorData(sites{iSite}(1));
        siteVertex(iSite,:) = det.V';

    ###########################################################################
    #    If only one site is present, use only one point for the error box.
    ###########################################################################
    if nSites == 1
        # ---- Use first input sky location.
        ra_ctr_deg = ra_ctr_deg(1);
        dec_ctr_deg = dec_ctr_deg(1);
        # ---- Command-line output.
        ra_search_deg = ra_ctr_deg;
        dec_search_deg = dec_ctr_deg;
        probabilities = 1;
        gridarea = []; #-- dummy value that can't be mistaken for physical number
        # ---- File output.
        # ---- Convert trigger ra,dec to phi,theta.
        [phi_ctr_rad, theta_ctr_rad] = radectoearth(ra_ctr_deg,dec_ctr_deg,gps);
        skyPositions = [theta_ctr_rad phi_ctr_rad 1 4*pi];
        dlmwrite(outputfile,skyPositions,'delimiter',' ','precision','#7.5f');     
        if verbose
            disp(['Single-site network: using only one grid point.']);
        end
        return
    end
    ###########################################################################
    #               Construct vector joining each pair of sites.
    ###########################################################################
    iPair = 1;
    iSite1 = 1;
    for iSite1 = 1:nSites-1
        for iSite2 = iSite1+1:nSites
           # ---- Construct structure listing names and displacements between
           #      pairs of sites.
           pairNames{iPair} = [sites{iSite1} sites{iSite2}];
           site1ToSite2(iPair,:) = siteVertex(iSite2,:) -  siteVertex(iSite1,:);
           # ---- Baseline in units of seconds.
           site1ToSite2_sec(iPair,:) = site1ToSite2(iPair,:) ./ speedOfLight; 
           iPair = iPair + 1;
        end
    end
    ###########################################################################
    #       Construct vector representing direction to error box centres
    #                          from Earth's centre.
    ###########################################################################
    # ---- Convert trigger ra,dec to phi,theta.
    [phi_ctr_rad, theta_ctr_rad] = radectoearth(ra_ctr_deg,dec_ctr_deg,gps);
    # ---- Construct the sky direction assumed from skyPosition_central.
    earthToGRB = [sin(theta_ctr_rad).*cos(phi_ctr_rad)... # x
                  sin(theta_ctr_rad).*sin(phi_ctr_rad)... # y
                  cos(theta_ctr_rad)];                    # z
    ###########################################################################
    #      For each pair of sites calculate angle with GRB line of sight.
    ###########################################################################
    # ---- For sites separated by distance delay_max [sec] have
    #        delay = delay_max cos(lambda)
    #      Differentiating yields
    #        ddelay = (-) alpha dlambda
    #      where
    #        alpha := delay_max sin(lambda).
    #      Here we find the baseline giving the largest value of alpha.
    #      This will produce the most conservative spacing of sky positions
    #      (dlambda) for a given tolerance in delay time error (ddelay).
    Npair = length(pairNames);
    for iPair = 1:Npair
        # ---- Some shorthand.
        v1 = site1ToSite2_sec(iPair,:).';
        v2 = earthToGRB;
       
        # ---- Check for 2 detectors at the same site (i.e., H1-H2).  In this
        #      case set alpha = 0; otherwise compute alpha.
        if (norm(v1)<1e-6)  #-- time delay of 1 microsec = 1000 foot baseline ...
           
            alpha(iPair) = 0;
            # ---- Dummy values for vebose output case.
            lambda_extrema = 0;
            lambda_opt = 0;
           
        else
            # ---- Invert dot product to calculate opening angle.
            lambda_ctr_rad = acos(v2*v1/norm(v1));
            # ---- Lambda should be between 0 and pi.
            if max(lambda_ctr_rad) > pi | min(lambda_ctr_rad) < 0
                disp(['lambda_ctr_rad: ']);
                    disp(lambda_ctr_rad);
                error('lambda_ctr_rad should be between 0 and pi.')
            end
            # ---- Add (subtract) skyPosError to find extreme angles.
            lambda_extrema  = [lambda_ctr_rad - skyPosError_rad, ...
                lambda_ctr_rad + skyPosError_rad];
            # ---- We are looking to maximise sin(abs(lambda)).
            # ---- If our range of lambda excludes pi/2 then choose the
            #      extrema closest to pi/2.
            [dummy,idx]=min(abs(pi/2-lambda_extrema),[],2);
            for ii = 1:length(idx)
                lambda_opt(ii) = lambda_extrema(ii,idx(ii));
            end
            # ---- If our range of lambda includes pi/2 this will maximise
            #      sin(abs(lambda)).
            ind = find(lambda_extrema(:,1) < pi/2 & lambda_extrema(:,2) > pi/2);
            if ~isempty(ind)
                lambda_opt(ind) = pi/2;
            end
            # ---- Compute alpha for this baseline.
            alpha(:,iPair) = norm(v1) * sin(lambda_opt(:));
           
        end
       
        # ---- Some output, if desired.
        if verbose
            fprintf(1,'For #s : \n', pairNames{iPair});
            fprintf(1,'lambda_min : #2.5f \n', min(lambda_extrema));
            fprintf(1,'lambda_max : #2.5f \n', max(lambda_extrema));
            fprintf(1,'lambda_opt : #2.5f \n', lambda_opt);
            fprintf(1,'alpha      : #2.5f \n', alpha(:,iPair));
        end
    end # -- Loop over pairs.
    clear v1 v2
    # ---- Find largest alpha value.
    alpha_max = max(alpha,[],2);
    # ---- Use alpha_max to calculate most conservative angularResolution which
    #      will be used to place out skyPositions.
    angularResolution = 2 * delay_tol ./ alpha_max;
    if verbose
        disp(['Desired angularResolution = ' num2str(angularResolution(:).') '.']);
    end
    ###########################################################################
    #                    Lay down grid of desired type.
    ###########################################################################
    switch lower(gridtype)
       
        case 'circular'
           
            # ---- Use separate circular grids for each error circle, removing
            #      overlap as needed.
            # ---- Number of error circles.
            Ncircle = size(earthToGRB,1);
            # ---- Number of points in each circular grid.
            Ngrid = zeros(Ncircle,1);
            # ---- Prepare storage.
            skyPositions = []; 
            # ---- Order in which to lay down circles: least dense to most
            #      dense.  (Gives smallest number of grid points, since we
            #      throw away points from later circles that overlap earlier
            #      circles.)  Reorder list of circles appropriately.
            [angularResolution,I] = sort(angularResolution,'descend');
            earthToGRB = earthToGRB(I,:);
            ra_ctr_deg = ra_ctr_deg(I);
            dec_ctr_deg = dec_ctr_deg(I);       
            skyPosError_rad = skyPosError_rad(I);
            skyPosError_deg = skyPosError_deg(I);
            sigma_deg = sigma_deg(I);
            # ---- Lay down grid for each error circle.
            for ii = 1:Ncircle
                # ---- Place skyPositions about North Pole in sky patch with
                #      radius = skyPosError_rad.
                #      3rd arg means we will place point at North Pole.
                #      4th arg means we will use 'optimal' spacing in theta.
                tempPositions = sinusoidalMap(angularResolution(ii),skyPosError_rad(ii),1,1);
                theta = tempPositions(:,1);
                phi = tempPositions(:,2);
                pOmega = tempPositions(:,3);
                dOmega = tempPositions(:,4);
                Ngrid(ii) = length(theta); 
                # ---- Rotate skyPositions to center them on each error circle.
                # ---- Convert from spherical coords to cartesian for convenience.
                V = CartesianPointingVector(phi,theta);
                # ---- The first vector in V should point to the North Pole.
                #      Figure out opening angle betwen this and the GRB.
                v1 = V(1,:);
                if dot(v1,[0,0,1])~=1
                    error('First point in non-rotated skyPositions map should be North Pole');
                end
               
                v2 = earthToGRB(ii,:);
                # ---- Construct vector which we will rotate skyPositions around.
                v3 = cross(v1,v2);
                # ---- Angle which we will rotate skyPositions about v3 by.
                psi = acos(dot(v1,v2)/(norm(v1)*norm(v2)));
                # ---- Loop over all skyPositions rotating them by psi about v3.
                Vrot = zeros(Ngrid(ii),3);
                for iVec = 1:size(V,1)
                    Vrot(iVec,:) = RotateVector(V(iVec,:),v3,psi);   
                end
                # ---- Convert rotated skyPositions back to spherical coordinates
                #      and replace original skyPositions with rotated ones.
                skyPositions = [skyPositions; ...
                    acos(Vrot(:,3)) atan2(Vrot(:,2),Vrot(:,1)) pOmega dOmega];
            end
           
            # ---- Compute angle between each grid point and all circle centers
            #      (needed for probability calculation).
            V = CartesianPointingVector(skyPositions(:,2),skyPositions(:,1));
            # ---- Location of error circle centers.
            [phi_ctr_rad, theta_ctr_rad] = radectoearth(ra_ctr_deg,dec_ctr_deg,gps);
            C = CartesianPointingVector(phi_ctr_rad,theta_ctr_rad);
            # ---- cos(Angle between grid points and error circle centers)
            overlap = V*C';
           
            # ---- Remove any grid points from circle N that lie within any of
            #      circles 1, ..., N-1, since they are redundant.
            for ii = Ncircle:-1:2
                index = cumsum(Ngrid);
                index = [index(ii-1)+1:index(ii)];
                drop = find(sum(overlap(index,1:ii-1)>repmat(cos(skyPosError_rad(1:ii-1).'),length(index),1),2));
                skyPositions(index(drop),:) = [];
                overlap(index(drop),:) = [];
            end
           
        case 'healpix'
           
            # ---- Use healpix to generate all-sky grid with desired
            #      angularResolution, then discard points outside error
            #      circles.  Gives a nice regular grid with meaningful area for
            #      each point, but density may be up to 4 x higher than needed.
            # ---- Number of points to cover sky at desired resolution -- based
            #      on finest angular resolution needed for all error boxes.
            Ntotal = 4*pi/min(angularResolution)^2;
            # ---- Make smallest healpix grid with this many points.
            n = [0:7];
            Nhealpix = 12*2.^(2*n);
            ind = find(Nhealpix>=Ntotal);
            if verbose
                disp(['Using healpix grid with n = ' num2str(n(ind(1))) ' (' ...
                    num2str(Nhealpix(ind(1))) ' sky positions, angular resolution = ' ...
                    num2str((4*pi/Nhealpix(ind(1))).^0.5) ').']);
            end
            [coordinates, solidAngles, probabilities] = healpix(n(ind(1)));
            skyPositions = [coordinates, probabilities, solidAngles];
            # ---- Unit three-vectors pointing to each grid point.
            V = CartesianPointingVector(skyPositions(:,2),skyPositions(:,1));
            # ---- Location of error circle centers.
            [phi_ctr_rad, theta_ctr_rad] = radectoearth(ra_ctr_deg,dec_ctr_deg,gps);
            C = CartesianPointingVector(phi_ctr_rad,theta_ctr_rad);
            # ---- Throw away grid points outside both error circles.
            overlap = V*C';
            # pass = max(overlap >= cos(skyPosError_rad),[],2);
            pass = max(overlap >= cos(skyPosError_rad+angularResolution),[],2); #-- probably overkill
            index = find(pass);
            skyPositions = skyPositions(index,:);
            overlap = overlap(index,:); #-- keep for probability calculations
            V = V(index,:); #-- keep for plots
    end
    # ---- Some output, if desired.
    if verbose
        fprintf(1,'We place #d grid points.\n', size(skyPositions,1));
    end
    ###########################################################################
    #                   Calculate probability of sky positions.
    ###########################################################################
    # -- If provided error is absurdly large use a uniform distribution
    if any(sigma_deg > 360)
      warning(['The provided 1-sigma error is greater than 360 degrees, using ' ...
               'a uniform penalty instead of a Fisher penalty']);
      skyPositions(:,3) = 1;
    else
    # ---- We assume that the probabilities of sky positions follow a Fisher
    #      distribution in the angle from the center of each error circle. Note
    #      that for each grid point we add the Fisher probabilities for each
    #      error circle.
    # ---- Compute the approximate kappa to get the 68# containment for the
    #      Fisher distribution.
    sigma = sigma_deg(:).' * (pi/180);
    kappa = (0.66 * sigma).^(-2);
    # ---- Verify that the approximation is not too bad.
    p1sigma = 0.68;
    containmentRadius = acos((kappa+log(1-p1sigma+p1sigma*exp(-2*kappa)))./kappa);
    if any(abs(containmentRadius-sigma) > 0.3 * sigma)
        error(['Discrepency between desired and generated containment ' ...
            'radius is greater than 30#'])
    end
    if any(abs(containmentRadius-sigma) > 0.1 * sigma)
        warning(['Discrepency between desired and generated containment ' ...
            'radius is greater than 10#'])
    end
    # ---- Compute probabilities, summing over all error circles.  Note that
    #      fisherpdf returns pdf for "polar" angle out from source, so to get
    #      2D probability distribution we need to divide by
    #      2*pi*sin(polar_angle).  Finally, we rescale so that the peak value
    #      is unity.
    kappa = repmat(kappa,size(skyPositions,1),1);
    skyPositions(:,3) = sum(exp(kappa.*(overlap-1)),2);
    skyPositions(:,3) = skyPositions(:,3) / max(skyPositions(:,3));
    # skyPositions(:,3) =  sum(kappa./(1-exp(-2*kappa)).*exp(kappa.*(overlap-1)),2)/...
    #     (kappa./(1-exp(-2*kappa)));
    end # -- choice between Fisher and uniform distribution
    if verbose
        # ---- Plot for debugging.
        figure; set(gca,'fontsize',16)
        plot3(V(:,1),V(:,2),V(:,3),'b.');
        hold on
        plot3(1.03*C(:,1),1.03*C(:,2),1.03*C(:,3),'m.','markersize',20);
        axis([-1.1 1.1 -1.1 1.1 -1.1 1.1])
        axis square
        grid on
        xlabel('x');
        ylabel('y');
        zlabel('z');
    end
    ###########################################################################
    #                         Convert back to ra,dec.
    ###########################################################################
    # ---- Generate output arguments.
    [ra_search_deg, dec_search_deg] = earthtoradec(skyPositions(:,2),skyPositions(:,1),gps);
    probabilities = skyPositions(:,3);
    gridarea = skyPositions(:,4);
    ###########################################################################
    #                         Write output file.
    ###########################################################################
    if ~strcmp(outputfile,'None')
        dlmwrite(outputfile,skyPositions,'delimiter',' ','precision','#7.5f');
    end
    # ---- Done.
    return ra_search_deg, dec_search_deg, probabilities, gridarea
    """
