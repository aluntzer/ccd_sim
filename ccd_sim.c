/**
 * @file   ccd_sim.c
 * @author Armin Luntzer (armin.luntzer@univie.ac.at),
 * @date   2022
 *
 * @copyright GPLv2
 * This program is free software; you can redistribute it and/or modify it
 * under the terms and conditions of the GNU General Public License,
 * version 2, as published by the Free Software Foundation.
 *
 * This program is distributed in the hope it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * @brief a fast CCD simulator extracted from the FEE simulator code

 */

#include <stddef.h>
#include <stdint.h>
#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>

#include <math.h>

#include <fitsio.h>


static uint16_t *CCD;



/**
 *  CCD characteristics for simulation
 */

/* the number of ccd pixels for a given row or column (assumed square) */
#define CCD_PIX_PER_AX			4510

/* the actually used CCD image section */
#define CCD_IMG_SEC_ROWS		3791
#define CCD_IMG_SEC_COLS		2255

/* a matching bin mode */
#define CCD_BIN_FRAME_6x6_ROWS		 639
#define CCD_BIN_FRAME_6x6_COLS		 384


/* pixel responsitivity in µV/electron */
#define CDD_RESP_uV_e		7.0
/* the gain equivalent, currently in discrepancy to value above*/
#define CCD_ADC_GAIN		40.

/* next few characteristics are just guesses based on a
 * typical CMOS ADC 5V input range and the responsitivity above
 * full well capacity in e-, */
#define CCD_N_FWC		714e3
/* dark signal, e-/pix/s, also CCD270-ish */
#define CCD_DARK		0.5
/* dark signal non-uniformity, as above */
#define CCD_DAR_NONUNI		0.05
/* readout noise e- rms */
#define CCD_NOISE		20.

/* number of electrons generated per 1 eV,
 * we will assume linear behaviour */
#define e_PER_eV		(55.0/200.0)
/* ccd thickness in µm, we us that for our cosmics */
#define CCD_THICKNESS_um	16.0
/* a ccd side length in mm */
#define CCD_SIDE_mm		81.8
/* ccd is partially covered for readout, this is the height
 * of the imaging area in mm */
#define CCD_IMG_HEIGHT_mm	68.24
/* the side dimensions of a pixel in µmn */
#define PIXEL_LEN_um		((81.8 * 1000.) / CCD_PIX_PER_AX)
/* we have 16 bit ADC res (allegedly) */
#define PIX_SATURATION		((1 << 16) - 1)


/* probability of charge transfer inefficiency occuring in a column */
#define CTI_PROB		0.1
/* bleed rate percentage if a CTI happens */
#define CTI_BLEED		0.1
/* probability of multi-pixel hits; set for testing; reasonable value: 0.002 */
#define MULTIPIX_HIT_PROB	0.5

/* number of pre-computed dark samples for faster simulation */
#define DARK_SAMPLES		128
/* simulate dark current (0/1) */
#define CFG_SIM_DARK		0
/* number of pre-computed readout noise samples for faster simulation */
#define RD_NOISE_SAMPLES	128


/**
 * particle simulation characteristics for CCD
 */

/* restrict the incident angle towards the CCD plane of the solar wind */
#define SOLAR_WIND_EL_ANGLE_MAX	(10. / 180. * M_PI)
/* restrict the azimuth angle in the CCD plane of the solar wind */
#define SOLAR_WIND_AZ_ANGLE_MAX	(30. / 180. * M_PI)

/* maximum random scattering angle (degrees) for particle deflection simulation;
 * the higher this value, the less the total probability of a deflected
 * particle trail occuring in the CCD image
 */
#define RUTHERFORD_SCATTER_ANGLE_MAX	90

/**
 * the particle energy loss follows a type of Landau distribution,
 * but appears to be relatively uniform regardless of the particle energy
 * 0.23 keV/µm of silicone appears to be a sensible value for our purposes, see
 * doi 10.1088/1748-0221/6/06/p06013
 *
 * given our ccd thickness, this would amount to a deposition on the order
 * of ~10k e- per pixel, which is only a few % of its FWC, which really
 * isn't all that much
 */
#define PARTICLE_ENERGY_LOSS_PER_UM_EV 230.


/**
  * solar activity selector, range 0-1
 */

#define SOLAR_ACT		1.0


/**
 *  typical environment in earth neighbourhood (do not change unless you know
 *  what you're doing)
 */

/* values values from "The SMILE Soft X-ray Imager (SXI) CCD design and
 * development" (Soman et al) */
/* nominal photon energy range */
#define SWCX_PHOT_EV_MIN	 200.0
#define SWCX_PHOT_EV_MAX	2000.0

/* some values from SMILE SXI CCD Testing and Calibration Event
 * Detection Methodology TN 1.2 (Soman et al)
 */
/* the following integrated event rates are in counts/CCD/s for the illuminated
 * section of a CCD;
 *
 * note that we could add sigmas for the particular rates, but for now we'll
 * just use these as the mean with a sigma of 1
 */
/* solar wind exchange x-ray event rate for low and high solar activity */
#define SWCX_CCD_RATE_MIN	5.134
#define SWCX_CCD_RATE_MAX	82.150
/* soft x-ray thermal galacitc background event rate */
#define SXRB_CCD_RATE		15.403
/* (solar?) particle background event rate */
#define PB_CCD_RATE_MIN		0.627
#define PB_CCD_RATE_MAX		1.255
/* mean point source rate over FOV */
#define PS_CCD_RATE		0.657
/* 4pi cosmic ray flux (in particles per CCD as well) */
#define COSMIC_FLUX		24.61


/* we assume solar particle energies of 1 keV to 100 keV
 * taken from
 * Long-Term Fluences of Energetic Particles in the Heliosphere (Mewaldt et al)
 *
 * NOTE: the values for (what I assume) refer to solar wind contribution
 * (PB_CCD_RATE_...) appear a little low, on the other hand, their
 * direction more or less correspond to the plane fo the CCDs, since
 * we're looking mostly in a direction which is perpendicular-ish to the Sun.
 */
#define SOLAR_PARTICLE_EV_MIN		1e03
#define SOLAR_PARTICLE_EV_MAX		1e05

/* we assume cosmic particle energies of 10 MeV to 10 GeV */
#define COSMIC_PARTICLE_EV_MIN		1e07
#define COSMIC_PARTICLE_EV_MAX		1e11
/* energy distribution per magnitute over min (i.e. x=0, 1, 2 ...) */
#define PARTICLE_DROPOFF	1.0		/* extra drop */
#define PARTICLE_RATE_DROP(x)	powf(10., -(x) * PARTICLE_DROPOFF + 1.)





static void save_fits(const char *name, uint16_t *buf, long rows, long cols)
{
	long naxes[3];
	long n_elem;
	int status = 0;
	fitsfile *ff;
	long fpixel[3] = {1, 1, 1};



	naxes[0] = cols;
	naxes[1] = rows;
	naxes[2] = 1;
	n_elem = naxes[0] * naxes[1] * naxes[2];

	if (fits_create_file(&ff, name, &status)) {
		printf("error %d\n", status);
		exit(-1);
	}

	if (fits_create_img(ff, USHORT_IMG, 3, naxes, &status)) {
		printf("error %d\n", status);
		exit(-1);
	}

	if (fits_write_pix(ff, TUSHORT, fpixel, n_elem, buf, &status)) {
		printf("error %d (%d)\n",status, __LINE__);
		exit(-1);
	}

	if (fits_close_file(ff, &status)) {
		printf("error %d (%d)\n",status, __LINE__);
		exit(-1);
	}
}


/**
 * @brief get a random number between 0 and 1 following a logarithmic
 * distribution
 * @note the base resolution is fixed to 1/1000
 */
static float sim_rand_log(void)
{
#define LOG_RAND_RES 1000
	float u;
	const float res = 1.0 / LOG_RAND_RES;

	u = 1.0 - res * (rand() % LOG_RAND_RES);

	return log(u) / log(res);
}



/**
 * @brief a box-muller gaussian-like distribution
 */

static float sim_rand_gauss(void)
{
	static int gen;

	static float u, v;


	gen = !gen;

	if (!gen)
		return sqrtf(- 2.0 * logf(u)) * cosf(2.0 * M_PI * v);

	u = (rand() + 1.0) / (RAND_MAX + 2.0);
	v =  rand()        / (RAND_MAX + 1.0);

	return sqrtf(-2.0 * logf(u)) * sinf(2.0 * M_PI * v);
}


static uint16_t ccd_sim_get_swcx_ray(void)
{
	float p;

	/* we assume the incident x-ray energy is uniformly distributed */
	p = fmodf(rand(), (SWCX_PHOT_EV_MAX + 1 - SWCX_PHOT_EV_MIN) * 1000.) * 0.001;
        p += SWCX_PHOT_EV_MIN;
	p *= e_PER_eV;
	p *= CDD_RESP_uV_e;	/* scale to voltage-equivalent */

	if (p > PIX_SATURATION)
		return PIX_SATURATION;
	else
		return (uint16_t) p;
}


static void ccd_sim_add_swcx(uint16_t *frame, uint16_t tint_ms)
{
	size_t n = CCD_IMG_SEC_ROWS * CCD_IMG_SEC_COLS;
	const float sigma = 1.0;
	float amp;
	float ray;
	float tint = (float) tint_ms / 1000.0;

	size_t i;
	size_t x, y;
	size_t pix;


	/* our event amplitude is the integration time for the configured
	 * event rates
	 */
	amp = tint * (SWCX_CCD_RATE_MIN + SOLAR_ACT * (SWCX_CCD_RATE_MAX - SWCX_CCD_RATE_MIN));

	/* use the square of the amplitude to scale the noise */
	amp = amp + sqrtf(amp) * sigma * sim_rand_gauss();

	printf("SWCX: %d rays produced\n", (int) amp/2);

	for (i = 0; i < ((size_t) amp / 2); i++) {

		ray = ccd_sim_get_swcx_ray();
		pix = rand() % (n + 1);

		/* trigger on random occurence (retval != 0) */
		if (rand() % ((int) (1. / MULTIPIX_HIT_PROB))) {
			frame[pix] += ray;
		} else {
			float fray = (float) ray;
			x = pix % CCD_IMG_SEC_COLS;
			y = (pix - x) / CCD_IMG_SEC_COLS;

			/* XXX meh... this function needs cleanup
			 * what's going on here: distribute hit power
			 * to adjacent pixels
			 */
			while (fray > 0.0) {
				int yy = 1 - rand() % 2;
				int xx = 1 - rand() % 2;
				ssize_t pp = (yy + y) * CCD_IMG_SEC_COLS + (xx + x);
				float bleedoff = ((float) (rand() % 100)) * 0.01 * fray;

				/* out of bounds? */
				if (pp < 0 || pp > (ssize_t) n)
					continue;

				/* make sure to reasonably abort the loop */
				if (bleedoff < 0.05 * (float) ray)
					bleedoff = fray;

				frame[pp] += bleedoff;
				fray -= bleedoff;
			}
		}
	}
}



static void ccd_sim_add_sxrb(uint16_t *frame, uint16_t tint_ms)
{
	size_t n = CCD_IMG_SEC_ROWS * CCD_IMG_SEC_COLS;
	const float sigma = 1.0;
	float amp;
	float tint = (float) tint_ms / 1000.0;

	size_t i;

	/* XXX add configurable CTI effect */

	/* our event amplitude is the integration time for the configured
	 * event rates
	 */
	amp = tint * SXRB_CCD_RATE;

	/* use the square of the amplitude to scale the noise */
	amp = amp + sqrtf(amp) * sigma * sim_rand_gauss();
	/* we use the same energies as SWCX */
	for (i = 0; i < ((size_t) amp / 2); i++)
		frame[rand() % (n + 1)] += ccd_sim_get_swcx_ray();
}


/**
 * @brief get a particle within the given energy range distribution
 */

static float ccd_sim_get_cosmic_particle(void)
{
	const float pmin = log10f(COSMIC_PARTICLE_EV_MIN);
	const float pmax = log10f(COSMIC_PARTICLE_EV_MAX);

	float r, p;

	/* get a energy range exponent */
	r = fmodf(rand(), (pmax + 1. - pmin) * 1000.) * 0.001;
	/* distribute to (logarithmic) particle rate */
	r = PARTICLE_RATE_DROP(r);

	/* get energy of the particle */
	p = COSMIC_PARTICLE_EV_MIN + (COSMIC_PARTICLE_EV_MAX - COSMIC_PARTICLE_EV_MIN) * r;
	p *= e_PER_eV;

	return p;
}


/**
 * @brief get a solar particle within the given energy range distribution
 */

static float ccd_sim_get_solar_particle(void)
{
	const float pmin = log10f(SOLAR_PARTICLE_EV_MIN);
	const float pmax = log10f(SOLAR_PARTICLE_EV_MAX);

	float r, p;

	/* get a energy range exponent */
	r = fmodf(rand(), (pmax + 1. - pmin) * 1000.) * 0.001;

	/* we assume equal probability for solar wind components */

	/* get energy of the particle */
	p = SOLAR_PARTICLE_EV_MIN + (SOLAR_PARTICLE_EV_MAX - SOLAR_PARTICLE_EV_MIN) * r;
	p *= e_PER_eV;

	return p;
}


/**
 * @brief get the scattering fraction of protons on Si atoms for
 *	  a CCD pixel
 *
 * @param p_eV the energy of the particle
 * @param theta the (minimum) scattering angle for the fraction scattered
 *
 * @note Rutherford scattering requires the target to only be a few µm thick
 *	 and preferably uses alpha particles as projectiles.
 *	 We're only interested in an approximate behaviour, so we can
 *	 easily get away with any deviations.
 */


static float ccd_sim_get_scatter_fraction(float p_eV, float theta)
{
	float t, r;
	float sigma;
	float f;

	/* we select between hydrogen and helium cores (~8%) */
	const float zp = (float) ((rand() % 100) <= 8 ? 2 : 1);
	const float Z  = 14.;				/* atomic number of Si */
	const float A  = 2.* Z;				/* mass number of  Si */
	const float rho= 2.33 * 1000.;			/* density of Si (kg/m^3) */
	const float Zp = powf(zp, 2.) / 4.;		/* projectile charge factor */
	const float k  = 8.9875517923e9;		/* Coloumb's constant */
	const float e  = -1.602176634e-19;		/* electron charge */
	const float eV = -e;				/* 1 eV in J */
	const float NA = 6.02214076e23;			/* Avogadro's number */
	const float L  = CCD_THICKNESS_um * 1e-6;	/* target thickness */


	t = k * pow(e, 2.) / (p_eV * eV);
	r = (1. + cosf(theta)) / (1. - cosf(theta));
	sigma = M_PI * Zp * pow(Z, 2.) * pow(t, 2.) * r;

	f = (NA * L * rho * sigma) / (A * 1e-3);

	if (f > 1.0)
		return 1.0;

	return f;
}

/**
 * @brief create particle traces in the CCD
 *
 * @param tint_ms the integration time in ms
 * @param solar 0 = cosmics, 1 = solar
 *
 */

static void ccd_sim_add_particles(uint16_t *frame, uint16_t tint_ms, int solar)
{
	size_t n = CCD_IMG_SEC_ROWS * CCD_IMG_SEC_COLS;
	const float sigma = 1.0;
	float amp;
	float tint = (float) tint_ms / 1000.0;

	size_t i;

	float x, y;
	float p_ev;
	float phi, theta;
	float d, r;
	float dx, dy;
	float d_ev;


	float *ccd;

	float deflection_angle;
	unsigned int deflection_rate;


	ccd = (float *) calloc(sizeof(float), n);
	if (!ccd) {
		perror("malloc");
		exit(-1);
	}

	/* our event amplitude is the integration time for the configured
	 * event rates
	 */
	if (solar)
		amp = tint * (PB_CCD_RATE_MIN + (PB_CCD_RATE_MAX - PB_CCD_RATE_MIN) * SOLAR_ACT);
	else
		amp = tint * COSMIC_FLUX;

	/* use the square of the amplitude to scale the noise */
	amp = amp + sqrtf(amp) * sigma * sim_rand_gauss();

	for (i = 0; i < ((size_t) amp / 2); i++) {


		/* initial particle energy */
		if (solar)
			p_ev = ccd_sim_get_solar_particle();
		else
			p_ev = ccd_sim_get_cosmic_particle();



		/* pixel of particle entry into CCD */
		x = (float) (rand() % (CCD_IMG_SEC_COLS + 1));
		y = (float) (rand() % (CCD_IMG_SEC_ROWS + 1));


		/* angle from CCD plane */
		if (solar) /* shallow */
			phi = fmod(rand(), SOLAR_WIND_EL_ANGLE_MAX * 1000.) * 0.001;
		else
			phi = fmod(rand(), M_PI_2 * 1000. ) * 0.001;

		/* direction within plane */
		if (solar) /* just one side */
			theta = 0.5 * SOLAR_WIND_AZ_ANGLE_MAX - fmod(rand(),  SOLAR_WIND_AZ_ANGLE_MAX * 1000.) * 0.001;
		else	/* anywhere */
			theta = M_PI - fmod(rand(), 2. * M_PI * 1000.) * 0.001;


restart:
		/* max distance traveled through ccd */
		d = CCD_THICKNESS_um / tanf(phi);

		/* get a random deflection angle */
#if 1
		/* log distribution */
		deflection_angle = sim_rand_log() * (RUTHERFORD_SCATTER_ANGLE_MAX / 180. * M_PI);
#else
		/* uniform */
		deflection_angle = ((float) (1 + rand() % RUTHERFORD_SCATTER_ANGLE_MAX )) / 180. * M_PI;
#endif
		/* our rate for rand() */
		deflection_rate =  (unsigned int) (1.0 / ccd_sim_get_scatter_fraction(p_ev, deflection_angle));
#if 0
		printf("deflection rate %u angle %g frac: %g ev %g\n", deflection_rate, deflection_angle / M_PI * 180., ccd_sim_get_scatter_fraction(p_ev, deflection_angle), p_ev);
#endif
		/* step size in x and y direction, we compute
		 * one pixel at a time and assume they are all cubes,
		 * so the max diagonal is sqrt(2) * thickness for the
		 * absorption
		 */
		dx = (PIXEL_LEN_um * M_SQRT2) * sinf(theta);
		dy = (PIXEL_LEN_um * M_SQRT2) * cosf(theta);
		/* max distance per pixel in vertical direction */
		r  = (CCD_THICKNESS_um * M_SQRT2) * cosf(phi);

		/* scale energy loss by longest distance travelelled
		 * within a pixel and the material-specific loss
		*/
		if (fabsf(dx) > fabsf(dy))
			d_ev = fabsf(dx);
		else
			d_ev = fabsf(dy);

		/* in case vertical travel component is the longest */
		if (fabsf(d_ev) < fabsf(r))
			d_ev = fabsf(r);

		/* scale to material */
		d_ev *= PARTICLE_ENERGY_LOSS_PER_UM_EV;

#if 0
		printf("energy: %g, x: %g y: %g, phi %g, theta %g streak %g, dx %g, dy %g, r %g, d_ev keV %g\n",
		       p_ev, x, y, phi * 180. / M_PI, theta * 180. / M_PI, d, dx, dy, r, d_ev/1000);
#endif

		/* XXX scale back to integer CCD pixels (not pretty..) */
		dx /= PIXEL_LEN_um * M_SQRT2;
		dy /= PIXEL_LEN_um * M_SQRT2;


		while (1) {

			size_t pix;

			if (x <= 0)
				break;
			if (y <= 0)
				break;
			if (x > CCD_IMG_SEC_COLS)
				break;
			if (y > CCD_IMG_SEC_ROWS)
				break;
			if (d < 0.)
				break;
			/* could set configurable cutoff here */
			if (p_ev < 0.)
				break;

			pix = (size_t) y * CCD_IMG_SEC_COLS + (size_t) x;
			ccd[pix] += d_ev * e_PER_eV * CDD_RESP_uV_e; // * (p_ev/p0) * logf(d/d0);

			x += dx;
			y += dy;

			p_ev -= d_ev;
			d -= r;

			if (rand() % (deflection_rate + 1) == 0) {
				float ratio = ((float) (rand(), 50)) * 0.01;
				float sign = (rand() & 0x1 ? -1.0 : 1.0);

				/* deflect the sucker */
				theta += sign * deflection_angle * ratio;
				sign = (rand() & 0x1 ? -1.0 : 1.0);
				phi   += sign * deflection_angle * (1. - ratio);
				goto restart;
			}

		}

	}

	for (i = 0; i < n; i++) {

		float tot = frame[i] + ccd[i];

		if (tot < (float) PIX_SATURATION)
			frame[i] = (uint16_t) tot;
		else /* else saturate, TODO: bleed charges (maybe) */
			frame[i] = PIX_SATURATION;
	}


	free(ccd);
}



/**
 * we don't really need a dark sim, the effective amplitude variation is way
 * to low to be significant
 */
__attribute__((unused))
static void ccd_sim_add_dark(uint16_t tint_ms)
{
	float *noisearr;
	size_t n = CCD_IMG_SEC_ROWS * CCD_IMG_SEC_COLS;
	float amp;
	float tint = (float) tint_ms / 1000.0;

	size_t i;



	/* total average accumulated dark current amplitude */
	amp = tint * CCD_DARK;

	/* use the square of the amplitude to scale the noise */
	amp = amp + sqrtf(amp);

	noisearr = (float *) calloc(sizeof(float), DARK_SAMPLES);
	if (!noisearr) {
		perror("malloc");
		exit(-1);
	}

	for (i = 0; i < DARK_SAMPLES; i++) {
		noisearr[i] =  amp + fmodf(sim_rand_gauss(), CCD_DAR_NONUNI * 1000.) * 0.001;
		noisearr[i]*= CDD_RESP_uV_e;	/* scale to voltage-equivalent */
	}


	/* fill CCD */
	for (i = 0; i < n; i++)
		CCD[i] += noisearr[rand() % (DARK_SAMPLES + 1)];

	free(noisearr);
}


static void ccd_sim_add_rd_noise(uint16_t *ccd, size_t n)
{
	float *noisearr;
	const float sigma = 1.0;
	float amp;

	size_t i;

	struct timeval t0, t;
	double elapsed_time;

	gettimeofday(&t0, NULL);
	/* total average accumulated dark current amplitude */
	amp = CCD_NOISE;

	noisearr = (float *) calloc(sizeof(float), RD_NOISE_SAMPLES);
	if (!noisearr) {
		perror("malloc");
		exit(-1);
	}

	for (i = 0; i < RD_NOISE_SAMPLES; i++) {
		/* use the square of the amplitude to scale the noise */
		noisearr[i] =  amp + sqrtf(amp) * sigma * sim_rand_gauss();
		noisearr[i]*= CDD_RESP_uV_e;	/* scale to voltage-equivalent */
	}


	/* add noise */
	for (i = 0; i < n; i++)
		ccd[i] += noisearr[rand() % (RD_NOISE_SAMPLES + 1)];

	free(noisearr);

	/* time elapsed in ms */
	gettimeofday(&t, NULL);
	elapsed_time  = (t.tv_sec  - t0.tv_sec)  * 1000.0;
	elapsed_time += (t.tv_usec - t0.tv_usec) / 1000.0;
	printf("readout noise in %g ms\n", elapsed_time);
}




static void ccd_sim_clear(void)
{
	size_t n = CCD_IMG_SEC_ROWS * CCD_IMG_SEC_COLS * sizeof(uint16_t);

	memset(CCD, 0, n);
}

static void ccd_sim_refresh(uint16_t tint_ms)
{
	struct timeval t0, t;
	double elapsed_time;


	gettimeofday(&t0, NULL);

	ccd_sim_clear();

	if (CFG_SIM_DARK)
		ccd_sim_add_dark(tint_ms);

	ccd_sim_add_swcx(CCD, tint_ms);
	ccd_sim_add_sxrb(CCD, tint_ms);

	/* solar and cosmic particles */
	ccd_sim_add_particles(CCD, tint_ms, 0);
	ccd_sim_add_particles(CCD, tint_ms, 1);

	/* time elapsed in ms */
	gettimeofday(&t, NULL);
	elapsed_time  = (t.tv_sec  - t0.tv_sec)  * 1000.0;
	elapsed_time += (t.tv_usec - t0.tv_usec) / 1000.0;
	printf("ccd refresh in %g ms\n", elapsed_time);

}



/**
 * @brief extract ccd data for frame transfer mode
 *
 * @param ccd  the ccd buffer
 * @param rows the rebinned rows in the output
 * @param cols the rebinned columns in the output
 * @param bins the pixel square that is being rebinned from the CCD
 *
 * @returns the allocated rebinned buffer
 */

static uint16_t *ccd_sim_get_bin_data(uint16_t *ccd, size_t rows, size_t cols,
				      size_t bins)
{
	size_t i, j;
	size_t x, y;
	size_t rw, cl;

	uint16_t *buf;
	uint16_t *acc;

	struct timeval t0, t;
	double elapsed_time;

	/* we need a cleared buffer to accumulate the bins */
	buf = calloc(sizeof(uint16_t), rows * cols);
	if (!buf) {
		perror("malloc");
		exit(-1);
	}

	if (bins == 1) {
		memcpy(buf, ccd, rows * cols * sizeof(uint16_t));
		return buf;
	}

	/* our accumulator */
	acc = malloc(CCD_IMG_SEC_COLS * sizeof(uint16_t));
	if (!acc) {
		perror("malloc");
		exit(-1);
	}

	gettimeofday(&t0, NULL);

	/* the binned data contain overscan in the real FEE, i.e. the "edge"
	 * pixels contain CCD bias values, but we just ignore those
	 * and round down to the next integer
	 * note that we keep the nominal size, we just don't fill the
	 * out-of-bounds samples with values
	 */
	rw = CCD_IMG_SEC_ROWS / bins;
	cl = CCD_IMG_SEC_COLS / bins;


	/* rebinned rows */
	for (y = 0; y < rw; y++) {

		/* clear line accumulator */
		memset(acc, 0x0, CCD_IMG_SEC_COLS * sizeof(uint16_t));

		/* accumulate the data lines in "blocks" of bins */
		for (i = 0; i < bins; i++) {

			/* the offset into the buffer of the start of the
			 * next row; note that the original number of rows
			 * is required, otherwise the data would be skewed
			 */
			size_t y0 = (y * bins + i) * CCD_IMG_SEC_COLS;

			for (j = 0; j < CCD_IMG_SEC_COLS; j++)
				acc[j] += ccd[y0 + j];

		}

		/* accumulate blocks of bins into columns */
		for (x = 0; x < cl; x++) {

			for (i = 0; i < bins; i++)
				buf[y * cols + x] += acc[x * bins + i];

		}

	}

	/* time in ms  */
	gettimeofday(&t, NULL);
	elapsed_time  = (t.tv_sec  - t0.tv_sec)  * 1000.0;
	elapsed_time += (t.tv_usec - t0.tv_usec) / 1000.0;
	printf("rebinned in %g ms\n", elapsed_time);

	free(acc);

	return buf;
}




static void ccd_sim_with_rebin(uint16_t tint_ms)
{
	/* these are taken from the original FEE mode */
	size_t rows = CCD_BIN_FRAME_6x6_ROWS;
	size_t cols = CCD_BIN_FRAME_6x6_COLS;
	size_t bins = 6;

	uint16_t *BIN;

	ccd_sim_refresh(tint_ms);

	BIN = ccd_sim_get_bin_data(CCD, rows, cols, bins);
	ccd_sim_add_rd_noise(BIN, rows * cols);


	/* add read noise to the full frame CCD */
	ccd_sim_add_rd_noise(CCD, CCD_IMG_SEC_ROWS * CCD_IMG_SEC_COLS);
	save_fits("!CCD.fits", CCD, CCD_IMG_SEC_ROWS, CCD_IMG_SEC_COLS);
	save_fits("!BIN.fits", BIN, rows, cols);


	free(BIN);
}



int main(void)
{
	const size_t img_side_pix = CCD_IMG_SEC_ROWS * CCD_IMG_SEC_COLS;

	CCD = (uint16_t *) malloc(img_side_pix * sizeof(uint16_t));
	if (!CCD) {
		perror("malloc");
		exit(-1);
	}

	/* 4s integration time */
	ccd_sim_with_rebin(4000);

	free(CCD);

	return 0;
}
