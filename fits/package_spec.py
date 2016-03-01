import numpy as np
import glob
import pyfits
import os
from lsst.sims.photUtils import Sed,Bandpass
from lsst.utils import getPackageDir
import healpy as hp

def packageAirglow():
    dataDir = getPackageDir('SIMS_SKYBRIGHTNESS_DATA')
    outDir = os.path.join(dataDir, 'ESO_Spectra/Airglow')

    # Read in all the spectra from ESO call and package into a single npz file
    files = glob.glob('Airglow/skytable*.fits')

    temp = pyfits.open(files[0])
    wave = temp[1].data['lam'].copy()*1e3

    airmasses = []
    solarFlux = []
    specs = []

    for i,filename in enumerate(files):
        fits = pyfits.open(filename)
        if np.max(fits[1].data['flux']) > 0:
            specs.append(fits[1].data['flux'].copy())
            header = fits[0].header['comment']
            for card in header:
                if 'SKYMODEL.TARGET.AIRMASS' in card:
                    airmasses.append(float(card.split('=')[-1]))
                elif 'SKYMODEL.MSOLFLUX' in card:
                    solarFlux.append(float(card.split('=')[-1]))

    airmasses = np.array(airmasses)
    solarFlux = np.array(solarFlux)

    nrec = airmasses.size
    nwave = wave.size

    dtype = [('airmass', 'float'),
             ('solarFlux', 'float'),
             ('spectra', 'float', (nwave)), ('mags', 'float', (6))]
    Spectra = np.zeros(nrec, dtype=dtype)
    Spectra['airmass'] = airmasses
    Spectra['solarFlux'] = solarFlux
    Spectra['spectra'] = specs

    hPlank = 6.626068e-27 # erg s
    cLight = 2.99792458e10 # cm/s

    # Convert spectra from ph/s/m2/micron/arcsec2 to erg/s/cm2/nm/arcsec2
    Spectra['spectra'] = Spectra['spectra']/(100.**2)*hPlank*cLight/(wave*1e-7)/1e3

    # Sort things since this might be helpful later
    Spectra.sort(order=['airmass','solarFlux'])

    # Load LSST filters
    throughPath = os.path.join(getPackageDir('throughputs'),'baseline')
    keys = ['u','g','r','i','z','y']
    nfilt = len(keys)
    filters = {}
    for filtername in keys:
        bp = np.loadtxt(os.path.join(throughPath, 'filter_'+filtername+'.dat'),
                        dtype=zip(['wave','trans'],[float]*2 ))
        tempB = Bandpass()
        tempB.setBandpass(bp['wave'],bp['trans'])
        filters[filtername] = tempB

    filterWave = np.array([filters[f].calcEffWavelen()[0] for f in keys ])

    for i,spectrum in enumerate(Spectra['spectra']):
        tempSed = Sed()
        tempSed.setSED(wave,flambda=spectrum)
        for j,filtName in enumerate(keys):
            try:
                Spectra['mags'][i][j] = tempSed.calcMag(filters[filtName])
            except:
                pass

    np.savez(os.path.join(outDir,'airglowSpectra.npz'), wave = wave, spec=Spectra, filterWave=filterWave)


def packageLowerAtm():

    dataDir = getPackageDir('SIMS_SKYBRIGHTNESS_DATA')
    outDir = os.path.join(dataDir, 'ESO_Spectra/LowerAtm')

    # Read in all the spectra from ESO call and package into a single npz file

    files = glob.glob('LowerAtm/skytable*.fits')

    temp = pyfits.open(files[0])
    wave = temp[1].data['lam'].copy()*1e3

    airmasses = []
    nightTimes = []
    specs = []

    for i,filename in enumerate(files):
        fits = pyfits.open(filename)
        if np.max(fits[1].data['flux']) > 0:
            specs.append(fits[1].data['flux'].copy())
            header = fits[0].header['comment']
            for card in header:
                if 'SKYMODEL.TARGET.AIRMASS' in card:
                    airmasses.append(float(card.split('=')[-1]))
                elif 'SKYMODEL.TIME' in card:
                    nightTimes.append(float(card.split('=')[-1]))

    airmasses = np.array(airmasses)
    nigtTimes = np.array(nightTimes)

    nrec = airmasses.size
    nwave = wave.size

    dtype = [('airmass', 'float'),
             ('nightTimes', 'float'),
             ('spectra', 'float', (nwave)), ('mags', 'float', (6))]
    Spectra = np.zeros(nrec, dtype=dtype)
    Spectra['airmass'] = airmasses
    Spectra['nightTimes'] = nightTimes
    Spectra['spectra'] = specs

    hPlank = 6.626068e-27 # erg s
    cLight = 2.99792458e10 # cm/s

    # Convert spectra from ph/s/m2/micron/arcsec2 to erg/s/cm2/nm/arcsec2
    Spectra['spectra'] = Spectra['spectra']/(100.**2)*hPlank*cLight/(wave*1e-7)/1e3

    # Sort things since this might be helpful later
    Spectra.sort(order=['airmass','nightTimes'])

    # Load LSST filters
    throughPath = os.path.join(getPackageDir('throughputs'),'baseline')
    keys = ['u','g','r','i','z','y']
    nfilt = len(keys)
    filters = {}
    for filtername in keys:
        bp = np.loadtxt(os.path.join(throughPath, 'filter_'+filtername+'.dat'),
                        dtype=zip(['wave','trans'],[float]*2 ))
        tempB = Bandpass()
        tempB.setBandpass(bp['wave'],bp['trans'])
        filters[filtername] = tempB

    filterWave = np.array([filters[f].calcEffWavelen()[0] for f in keys ])

    for i,spectrum in enumerate(Spectra['spectra']):
        tempSed = Sed()
        tempSed.setSED(wave,flambda=spectrum)
        for j,filtName in enumerate(keys):
            try:
                Spectra['mags'][i][j] = tempSed.calcMag(filters[filtName])
            except:
                pass

    np.savez(os.path.join(outDir,'Spectra.npz'), wave = wave, spec=Spectra, filterWave=filterWave)


def packageMoon():
    dataDir = getPackageDir('SIMS_SKYBRIGHTNESS_DATA')
    outDir = os.path.join(dataDir, 'ESO_Spectra/Moon')

    # Read in all the spectra from ESO call and package into a single npz file

    files = glob.glob('Moon/skytable*.fits')

    temp = pyfits.open(files[0])
    moonWave = temp[1].data['lam'].copy()*1e3

    moonSpec = []
    moonAM = [] # Actually target airmass
    moonAlt=[]  # altitude of the moon
    moonSunSep=[] # moon Phase
    moonTargetSep=[]
    hpid=[]
    for i,filename in enumerate(files):
        fits = pyfits.open(filename)
        try:
            if np.max(fits[1].data['flux']) > 0:
                moonSpec.append(fits[1].data['flux'].copy())
                header = fits[0].header['comment']
                for card in header:
                    if 'SKYMODEL.MOON.SUN.SEP' in card:
                        moonSunSep.append(float(card.split('=')[-1]))
                    elif 'SKYMODEL.TARGET.AIRMASS' in card:
                        #moonAM.append( 1./np.cos(np.radians(90.-float(card.split('=')[-1]))) )
                        moonAM.append( float(card.split('=')[-1]) )
                    elif 'SKYMODEL.MOON.TARGET.SEP' in card:
                        moonTargetSep.append(float(card.split('=')[-1]))
                    elif 'SKYMODEL.MOON.ALT' in card:
                        moonAlt.append(float(card.split('=')[-1]))
        except:
            print filename, ' Failed'

    import healpy as hp
    from lsst.sims.utils import haversine

    nside = 4
    lat, az = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    alt = np.pi/2.-lat
    airmass = 1./np.cos(np.pi/2.-alt)

    # Only need low airmass and then 1/2 to sky
    good = np.where( (az >= 0) & (az <= np.pi) & (airmass <=2.6) & (airmass >= 1.) )
    airmass = airmass[good]
    alt=alt[good]
    az = az[good]

    moonAM = np.array(moonAM)
    moonAlt = np.array(moonAlt)
    moonSunSep = np.array(moonSunSep)
    moonTargetSep = np.array(moonTargetSep)
    moonAzDiff = moonTargetSep*0
    targetAlt = np.pi/2.-np.arccos(1./moonAM)
    # Compute the azimuth difference given the moon-target-seperation
    # Let's just do a stupid loop:
    for i in np.arange(targetAlt.size):
        possibleDistances = haversine(0., np.radians(moonAlt[i]),  az, az*0+targetAlt[i])
        diff = np.abs(possibleDistances - np.radians(moonTargetSep[i]))
        good = np.where(diff == diff.min())
        moonAzDiff[i] = az[good][0]
        # ok, now I have an alt and az, I can convert that back to a healpix id.

        hpid.append(hp.ang2pix(nside, np.pi/2.-targetAlt[i], moonAzDiff[i]))

    nrec = moonAM.size
    nwave = moonWave.size

    dtype = [('hpid', 'int'),
             ('moonAltitude', 'float'),
             ('moonSunSep', 'float'),
             ('spectra', 'float', (nwave)), ('mags', 'float', (6))]
    moonSpectra = np.zeros(nrec, dtype=dtype)
    moonSpectra['hpid'] = hpid
    moonSpectra['moonAltitude'] = moonAlt
    moonSpectra['moonSunSep'] = moonSunSep
    moonSpectra['spectra'] = moonSpec

    hPlank = 6.626068e-27  # erg s
    cLight = 2.99792458e10 # cm/s

    # Convert spectra from ph/s/m2/micron/arcsec2 to erg/s/cm2/nm/arcsec2
    moonSpectra['spectra'] = moonSpectra['spectra']/(100.**2)*hPlank*cLight/(moonWave*1e-7)/1e3

    # Sort things since this might be helpful later
    moonSpectra.sort(order=['moonSunSep','moonAltitude', 'hpid'])

    # Crop off the incomplete ones

    good =np.where((moonSpectra['moonAltitude'] >= 0) & (moonSpectra['moonAltitude'] < 89.) )
    moonSpectra = moonSpectra[good]

    # Load LSST filters
    throughPath = os.path.join(getPackageDir('throughputs'),'baseline')
    keys = ['u','g','r','i','z','y']
    nfilt = len(keys)
    filters = {}
    for filtername in keys:
        bp = np.loadtxt(os.path.join(throughPath, 'filter_'+filtername+'.dat'),
                        dtype=zip(['wave','trans'],[float]*2 ))
        tempB = Bandpass()
        tempB.setBandpass(bp['wave'],bp['trans'])
        filters[filtername] = tempB

    filterWave = np.array([filters[f].calcEffWavelen()[0] for f in keys ])

    for i,spectrum in enumerate(moonSpectra['spectra']):
        tempSed = Sed()
        tempSed.setSED(moonWave,flambda=spectrum)
        for j,filtName in enumerate(keys):
            try:
                moonSpectra['mags'][i][j] = tempSed.calcMag(filters[filtName])
            except:
                pass

    nbreak=5
    nrec = np.size(moonSpectra)
    for i in np.arange(nbreak):
        np.savez(os.path.join(outDir,'moonSpectra_'+str(i)+'.npz'), wave = moonWave, spec=moonSpectra[i*nrec/nbreak:(i+1)*nrec/nbreak], filterWave=filterWave)


def packageScattered():

    dataDir = getPackageDir('SIMS_SKYBRIGHTNESS_DATA')
    outDir = os.path.join(dataDir, 'ESO_Spectra/ScatteredStarLight')

    # Read in all the spectra from ESO call and package into a single npz file

    files = glob.glob('ScatteredStarLight/skytable*.fits')

    temp = pyfits.open(files[0])
    wave = temp[1].data['lam'].copy()*1e3

    airmasses = []
    nightTimes = []
    specs = []

    for i,filename in enumerate(files):
        fits = pyfits.open(filename)
        if np.max(fits[1].data['flux']) > 0:
            specs.append(fits[1].data['flux'].copy())
            header = fits[0].header['comment']
            for card in header:
                if 'SKYMODEL.TARGET.AIRMASS' in card:
                    airmasses.append(float(card.split('=')[-1]))
                elif 'SKYMODEL.TIME' in card:
                    nightTimes.append(float(card.split('=')[-1]))


    airmasses = np.array(airmasses)
    nigtTimes = np.array(nightTimes)

    nrec = airmasses.size
    nwave = wave.size

    dtype = [('airmass', 'float'),
             ('nightTimes', 'float'),
             ('spectra', 'float', (nwave)), ('mags', 'float', (6))]
    Spectra = np.zeros(nrec, dtype=dtype)
    Spectra['airmass'] = airmasses
    Spectra['nightTimes'] = nightTimes
    Spectra['spectra'] = specs


    hPlank = 6.626068e-27 # erg s
    cLight = 2.99792458e10 # cm/s

    # Convert spectra from ph/s/m2/micron/arcsec2 to erg/s/cm2/nm/arcsec2
    Spectra['spectra'] = Spectra['spectra']/(100.**2)*hPlank*cLight/(wave*1e-7)/1e3

    # Sort things since this might be helpful later
    Spectra.sort(order=['airmass','nightTimes'])

    # Load LSST filters
    throughPath = os.path.join(getPackageDir('throughputs'),'baseline')
    keys = ['u','g','r','i','z','y']
    nfilt = len(keys)
    filters = {}
    for filtername in keys:
        bp = np.loadtxt(os.path.join(throughPath, 'filter_'+filtername+'.dat'),
                        dtype=zip(['wave','trans'],[float]*2 ))
        tempB = Bandpass()
        tempB.setBandpass(bp['wave'],bp['trans'])
        filters[filtername] = tempB

    filterWave = np.array([filters[f].calcEffWavelen()[0] for f in keys ])

    for i,spectrum in enumerate(Spectra['spectra']):
        tempSed = Sed()
        tempSed.setSED(wave,flambda=spectrum)
        for j,filtName in enumerate(keys):
            try:
                Spectra['mags'][i][j] = tempSed.calcMag(filters[filtName])
            except:
                pass

    np.savez(os.path.join(outDir,'scatteredStarLight.npz'), wave = wave, spec=Spectra, filterWave=filterWave)



def packageUpper():

    dataDir = getPackageDir('SIMS_SKYBRIGHTNESS_DATA')
    outDir = os.path.join(dataDir, 'ESO_Spectra/UpperAtm')

    # Read in all the spectra from ESO call and package into a single npz file

    files = glob.glob('UpperAtm/skytable*.fits')

    temp = pyfits.open(files[0])
    wave = temp[1].data['lam'].copy()*1e3

    airmasses = []
    nightTimes = []
    specs = []

    for i,filename in enumerate(files):
        fits = pyfits.open(filename)
        if np.max(fits[1].data['flux']) > 0:
            specs.append(fits[1].data['flux'].copy())
            header = fits[0].header['comment']
            for card in header:
                if 'SKYMODEL.TARGET.AIRMASS' in card:
                    airmasses.append(float(card.split('=')[-1]))
                elif 'SKYMODEL.TIME' in card:
                    nightTimes.append(float(card.split('=')[-1]))


    airmasses = np.array(airmasses)
    nigtTimes = np.array(nightTimes)

    nrec = airmasses.size
    nwave = wave.size

    dtype = [('airmass', 'float'),
             ('nightTimes', 'float'),
             ('spectra', 'float', (nwave)), ('mags', 'float', (6))]
    Spectra = np.zeros(nrec, dtype=dtype)
    Spectra['airmass'] = airmasses
    Spectra['nightTimes'] = nightTimes
    Spectra['spectra'] = specs


    hPlank = 6.626068e-27 # erg s
    cLight = 2.99792458e10 # cm/s

    # Convert spectra from ph/s/m2/micron/arcsec2 to erg/s/cm2/nm/arcsec2
    Spectra['spectra'] = Spectra['spectra']/(100.**2)*hPlank*cLight/(wave*1e-7)/1e3

    # Sort things since this might be helpful later
    Spectra.sort(order=['airmass','nightTimes'])

    # Load LSST filters
    throughPath = os.path.join(getPackageDir('throughputs'),'baseline')
    keys = ['u','g','r','i','z','y']
    nfilt = len(keys)
    filters = {}
    for filtername in keys:
        bp = np.loadtxt(os.path.join(throughPath, 'filter_'+filtername+'.dat'),
                        dtype=zip(['wave','trans'],[float]*2 ))
        tempB = Bandpass()
        tempB.setBandpass(bp['wave'],bp['trans'])
        filters[filtername] = tempB

    filterWave = np.array([filters[f].calcEffWavelen()[0] for f in keys ])

    for i,spectrum in enumerate(Spectra['spectra']):
        tempSed = Sed()
        tempSed.setSED(wave,flambda=spectrum)
        for j,filtName in enumerate(keys):
            try:
                Spectra['mags'][i][j] = tempSed.calcMag(filters[filtName])
            except:
                pass


    np.savez(os.path.join(outDir,'Spectra.npz'), wave = wave, spec=Spectra, filterWave=filterWave)


def packageZodiacal():

    dataDir = getPackageDir('SIMS_SKYBRIGHTNESS_DATA')
    outDir = os.path.join(dataDir, 'ESO_Spectra/Zodiacal')

    nside = 4

    # Read in all the spectra from ESO call and package into a single npz file

    files = glob.glob('Zodiacal/skytable*.fits')

    temp = pyfits.open(files[0])
    wave = temp[1].data['lam'].copy()*1e3

    airmasses = []
    eclLon = []
    eclLat = []
    specs = []

    for i,filename in enumerate(files):
        fits = pyfits.open(filename)
        if np.max(fits[1].data['flux']) > 0:
            specs.append(fits[1].data['flux'].copy())
            header = fits[0].header['comment']
            for card in header:
                if 'SKYMODEL.TARGET.AIRMASS' in card:
                    airmasses.append(float(card.split('=')[-1]))
                elif 'SKYMODEL.ECL.LON' in card:
                    eclLon.append(float(card.split('=')[-1]))
                elif 'SKYMODEL.ECL.LAT' in card:
                    eclLat.append(float(card.split('=')[-1]))


    airmasses = np.array(airmasses)
    eclLon = np.array(eclLon)
    eclLat = np.array(eclLat)

    wrapA = np.where(eclLon < 0.)
    eclLon[wrapA] = eclLon[wrapA]+360.

    uAM = np.unique(airmasses)
    nAM = uAM.size
    nwave = wave.size

    dtype = [('airmass', 'float'),
             ('hpid', 'int' ),
             ('spectra', 'float', (nwave)), ('mags', 'float', (6))]
    npix = hp.nside2npix(nside)
    Spectra = np.zeros(nAM*npix, dtype=dtype)
    for i,am in enumerate(uAM):
        Spectra['airmass'][i*npix:i*npix+npix] = am
        Spectra['hpid'][i*npix:i*npix+npix] = np.arange(npix)


    for am, lat, lon, spec in zip(airmasses,eclLat, eclLon, specs):
        hpid = hp.ang2pix(nside, np.radians(lat+90.), np.radians(lon) )
        good = np.where( (Spectra['airmass'] == am) & (Spectra['hpid'] == hpid))
        Spectra['spectra'][good] = spec.copy()

    hPlank = 6.626068e-27 # erg s
    cLight = 2.99792458e10 # cm/s

    # Convert spectra from ph/s/m2/micron/arcsec2 to erg/s/cm2/nm/arcsec2
    Spectra['spectra'] = Spectra['spectra']/(100.**2)*hPlank*cLight/(wave*1e-7)/1e3

    # Sort things since this might be helpful later
    Spectra.sort(order=['airmass', 'hpid'])

    # Load LSST filters
    throughPath = os.path.join(getPackageDir('throughputs'),'baseline')
    keys = ['u','g','r','i','z','y']
    nfilt = len(keys)
    filters = {}
    for filtername in keys:
        bp = np.loadtxt(os.path.join(throughPath, 'filter_'+filtername+'.dat'),
                        dtype=zip(['wave','trans'],[float]*2 ))
        tempB = Bandpass()
        tempB.setBandpass(bp['wave'],bp['trans'])
        filters[filtername] = tempB

    filterWave = np.array([filters[f].calcEffWavelen()[0] for f in keys ])

    for i,spectrum in enumerate(Spectra['spectra']):
        tempSed = Sed()
        tempSed.setSED(wave,flambda=spectrum)
        for j,filtName in enumerate(keys):
            try:
                Spectra['mags'][i][j] = tempSed.calcMag(filters[filtName])
            except:
                pass


    #span this over multiple files to store in github
    nbreak = 3
    nrec = np.size(Spectra)

    for i in np.arange(nbreak):
        np.savez(os.path.join(outDir,'zodiacalSpectra_'+str(i)+'.npz'), wave = wave,
                 spec=Spectra[i*nrec/nbreak:(i+1)*nrec/nbreak], filterWave=filterWave)


def packageMergedSpec():
    dataDir = getPackageDir('SIMS_SKYBRIGHTNESS_DATA')
    outDir = os.path.join(dataDir, 'ESO_Spectra/MergedSpec')

    # A large number of the background components only depend on Airmass, so we can merge those together
    npzs = ['LowerAtm/Spectra.npz',
            'ScatteredStarLight/scatteredStarLight.npz',
            'UpperAtm/Spectra.npz']
    files = [os.path.join(dataDir, 'ESO_Spectra', npz) for npz in npzs]


    temp = np.load(files[0])

    wave = temp['wave'].copy()
    spec = temp['spec'].copy()
    spec['spectra'] = spec['spectra']*0.
    spec['mags'] = spec['mags']*0.

    for filename in files:
        restored = np.load(filename)
        spec['spectra'] += restored['spec']['spectra']
        try:
            flux = 10.**(-0.4*(restored['spec']['mags']-np.log10(3631.)))
        except:
            import pdb ; pdb.set_trace()
        flux[np.where(restored['spec']['mags'] == 0.)] = 0.
        spec['mags'] += flux

    spec['mags'] = -2.5*np.log10(spec['mags'])+np.log10(3631.)

    np.savez(os.path.join(outDir,'mergedSpec.npz'), spec=spec, wave=wave, filterWave=temp['filterWave'])


if __name__ == '__main__':

    print 'packaging airglow'
    packageAirglow()
    print 'packaging lower atm'
    packageLowerAtm()
    print 'packaging moon'
    packageMoon()
    print 'packaging upper atm'
    packageUpper()
    print 'packaging scattered starlight'
    packageScattered()
    print 'packaging zodiacal'
    packageZodiacal()
    print 'packaging merged spectra'
    packageMergedSpec()
