import sys
import os
from socket import gethostname 
sys.path.append(os.curdir);
os.putenv("TDKPCONFPATH",os.curdir+"/tdkp_conf");
from negf import *
	
try:	
	cvar.fnames.init('create_structure');
	
	# extend path for module searching by current directory.
	# this is necessary when the given script file is a symbolic link
	sys.path.append(cvar.fnames.get_path());
	# TDKP configuration
	os.putenv("TDKPCONFPATH",cvar.fnames.get_path()+"conf");
	
	cvar.logmsg.add_device("std::cout");
	if cvar.mpi.get_rank()==constants.mpi_master_rank:
		cvar.logmsg.add_device(cvar.fnames.get_logfile());
	
	layers = [];
	# =====================================================================
	# START OF USER INPUT
	# =====================================================================
	PML = False;
	DATTADISS = False;
	RTD       = False;
	RTD_NELLO = False;
	QW        = False;
	PN        = False;
	IWCE      = False;
	NICHIA    = True;
	if DATTADISS:
		L  = 60;  # default: 96nm
		C  = 12;  # default: 12nm
		dx = 0.8; # default: 1.2nm
		mat = 'GaAs';
		#mat = 'Silicon';
		#structurename = 'structures/dattadiss'+str(L)+'nm'+str(dx)+'dx'+str(C);
		structurename = 'structures/dattadiss'+str(L)+'nm'+str(dx)+'dx';
		ndop  = 2e18;
		nplus = 2e19;
		#ndop  = 5e19;
		#nplus = 1e20;
		
		# syntax: name, material, xmole, length[nm], dx[nm], doping [cm-3]
		layers.append(['contact_0',mat, 0.0, 2*dx, dx, nplus]);
		layers.append(['left_nn'  ,mat, 0.0,    C, dx, nplus]);
		layers.append(['mid_n'    ,mat, 0.0,    L, dx, ndop]);
		layers.append(['right_nn' ,mat, 0.0,    C, dx, nplus]);
		layers.append(['contact_1',mat, 0.0, 2*dx, dx, nplus]);
	if RTD:
		structurename = 'structures/rtd';
		if PML:
			structurename = structurename+'_pml';
		
		# syntax: name, material, xmole, length[nm], dx[nm], doping [cm-3]
		ndop  = 2e18;
		nplus = 5e18;
		if PML:
			layers.append(['PML_mx'   ,   'GaAs', 0.0, 15.0, 0.2, ndop]);
			#layers.append(['left_nn'  ,   'GaAs', 0.0, 10.0, 1.0, nplus]);
			layers.append(['left_n'   ,   'GaAs', 0.0, 22.0, 0.2, ndop]);
			layers.append(['left_i'   ,   'GaAs', 0.0,  3.0, 0.2, 1e10]);
			layers.append(['barrier1' , 'AlGaAs', 0.3,  3.0, 0.2, 1e10]); 
			layers.append(['mid_bound',   'GaAs', 0.0,  4.0, 0.2, 1e10]);
			layers.append(['barrier2' , 'AlGaAs', 0.3,  3.0, 0.2, 1e10]);
			layers.append(['right_i'  ,   'GaAs', 0.0,  3.0, 0.2, 1e10]);
			layers.append(['right_n'  ,   'GaAs', 0.0, 22.0, 0.2, ndop]);
			#layers.append(['right_nn'  ,   'GaAs', 0.0, 10.0, 1.0, nplus]);
			layers.append(['PML_px'   ,   'GaAs', 0.0, 15.0, 0.2, ndop]);
		else:
			layers.append(['contact_0',   'GaAs', 0.0,  2.0, 1.0, ndop]);
			#layers.append(['left_nn'  ,   'GaAs', 0.0, 10.0, 1.0, nplus]);
			layers.append(['left_n'   ,   'GaAs', 0.0, 22.0, 1.0, ndop]);
			layers.append(['left_i'   ,   'GaAs', 0.0,  3.0, 1.0, 1e10]);
			layers.append(['barrier1' , 'AlGaAs', 0.3,  3.0, 0.5, 1e10]); 
			layers.append(['mid_bound',   'GaAs', 0.0,  4.0, 0.5, 1e10]);
			layers.append(['barrier2' , 'AlGaAs', 0.3,  3.0, 0.5, 1e10]);
			layers.append(['right_i'  ,   'GaAs', 0.0,  3.0, 1.0, 1e10]);
			layers.append(['right_n'  ,   'GaAs', 0.0, 22.0, 1.0, ndop]);
			#layers.append(['right_nn'  ,   'GaAs', 0.0, 10.0, 1.0, nplus]);
			layers.append(['contact_1',   'GaAs', 0.0,  2.0, 1.0, ndop]);
	if QW:
		structurename = 'structures/qw';
		if PML:
			structurename = structurename+'_pml';
		
		# syntax: name, material, xmole, length[nm], dx[nm], doping [cm-3]
		ndop  =  1e19;
		pdop  = -1e19;
		if PML:
			# fine discretization
#			layers.append(['PML_mx'   ,   'GaAs', 0.0, 15.0, 0.2,  pdop]);
#			layers.append(['left_p'   ,   'GaAs', 0.0, 10.0, 0.2,  pdop]);
#			layers.append(['left_i'   ,   'GaAs', 0.0,  7.0, 0.2,  1e10]);
#			layers.append(['barrier1' ,   'GaAs', 0.0,  3.0, 0.2,  1e10]); 
#			layers.append(['mid_bound', 'InGaAs', 0.3,  5.0, 0.2,  1e10]);
#			layers.append(['barrier2' ,   'GaAs', 0.0,  3.0, 0.2,  1e10]);
#			layers.append(['right_i'  ,   'GaAs', 0.0,  7.0, 0.2,  1e10]);
#			layers.append(['right_n'  ,   'GaAs', 0.0, 10.0, 0.2,  ndop]);
#			layers.append(['PML_px'   ,   'GaAs', 0.0, 15.0, 0.2,  ndop]);

			# same discretization as NEGF grid
			layers.append(['PML_mx'   ,   'GaAs', 0.0, 15.0, 0.2,  pdop]);
			layers.append(['left_p'   ,   'GaAs', 0.0, 10.0, 0.5,  pdop]);
			layers.append(['left_i'   ,   'GaAs', 0.0,  7.0, 0.5,  1e10]);
			layers.append(['barrier1' ,   'GaAs', 0.0,  3.0, 0.5,  1e10]); 
			layers.append(['iface1a'  ,   'GaAs', 0.0,  0.2, 0.1,  1e10]); 
			layers.append(['iface1b'  , 'InGaAs', 0.3,  0.2, 0.1,  1e10]); 
			layers.append(['mid_bound', 'InGaAs', 0.3,  5.0, 0.5,  1e10]); 
			layers.append(['iface2b'  , 'InGaAs', 0.3,  0.2, 0.1,  1e10]);
			layers.append(['iface2a'  ,   'GaAs', 0.0,  0.2, 0.1,  1e10]); 
			layers.append(['barrier2' ,   'GaAs', 0.0,  3.0, 1.0,  1e10]);
			layers.append(['right_i'  ,   'GaAs', 0.0,  7.0, 1.0,  1e10]);
			layers.append(['right_n'  ,   'GaAs', 0.0, 10.0, 1.0,  ndop]);
			layers.append(['PML_px'   ,   'GaAs', 0.0, 15.0, 0.2,  ndop]);
		else:
			layers.append(['contact_0',   'GaAs', 0.0,  1.0, 0.5,  pdop]);
			layers.append(['left_p'   ,   'GaAs', 0.0, 10.0, 0.5,  pdop]);
			layers.append(['left_i'   ,   'GaAs', 0.0,  7.0, 0.5,  1e10]);
			layers.append(['barrier1' ,   'GaAs', 0.0,  3.0, 0.5,  1e10]); 
			layers.append(['iface1a'  ,   'GaAs', 0.0,  0.2, 0.1,  1e10]); 
			layers.append(['iface1b'  , 'InGaAs', 0.3,  0.2, 0.1,  1e10]); 
			layers.append(['mid_bound', 'InGaAs', 0.3,  5.0, 0.5,  1e10]); 
			layers.append(['iface2b'  , 'InGaAs', 0.3,  0.2, 0.1,  1e10]);
			layers.append(['iface2a'  ,   'GaAs', 0.0,  0.2, 0.1,  1e10]); 
			layers.append(['barrier2' ,   'GaAs', 0.0,  3.0, 1.0,  1e10]);
			layers.append(['right_i'  ,   'GaAs', 0.0,  7.0, 1.0,  1e10]);
			layers.append(['right_n'  ,   'GaAs', 0.0, 10.0, 1.0,  ndop]);
			layers.append(['contact_1',   'GaAs', 0.0,  2.0, 1.0,  ndop]);
	if RTD_NELLO:
		structurename = 'structures/rtd_nello';
		if PML:
			structurename = structurename+'_pml';
		
		# syntax: name, material, xmole, length[nm], dx[nm], doping [cm-3]
		ndop  = 1e19;
		if PML:
			layers.append(['PML_mx'   ,   'GaAs', 0.0, 15.0, 0.2,  ndop]);
			layers.append(['left_n'   ,   'GaAs', 0.0, 40.0, 0.4,  ndop]);
			layers.append(['left_n2'  ,   'GaAs', 0.0,  3.0, 0.2,  ndop]);
			layers.append(['barrier1' , 'AlGaAs', 0.4,  3.0, 0.2,  1e10]); 
			layers.append(['mid_bound',   'GaAs', 0.0,  8.0, 0.2,  1e10]);
			layers.append(['barrier2' , 'AlGaAs', 0.4,  3.0, 0.2,  1e10]);
			layers.append(['right_n2' ,   'GaAs', 0.0,  3.0, 0.2,  ndop]);
			layers.append(['right_n'  ,   'GaAs', 0.0, 40.0, 0.4,  ndop]);
			layers.append(['PML_px'   ,   'GaAs', 0.0, 15.0, 0.2,  ndop]);
		else:
			layers.append(['contact_0',   'GaAs', 0.0,  2.5, 1.25, ndop]);
			layers.append(['left_n'   ,   'GaAs', 0.0, 40.0, 1.25, ndop]);
			layers.append(['left_n2'  ,   'GaAs', 0.0,  3.0,  0.5, ndop]);
			layers.append(['barrier1' , 'AlGaAs', 0.4,  3.0,  0.5, 1e10]); 
			layers.append(['mid_bound',   'GaAs', 0.0,  8.0,  0.5, 1e10]);
			layers.append(['barrier2' , 'AlGaAs', 0.4,  3.0,  0.5, 1e10]);
			layers.append(['right_n2' ,   'GaAs', 0.0,  3.0,  0.5, ndop]);
			layers.append(['right_n'  ,   'GaAs', 0.0, 40.0, 1.25, ndop]);
			layers.append(['contact_1',   'GaAs', 0.0,  2.5, 1.25, ndop]);
	if PN:
		# structurename = 'structures/pn';
		# structurename = 'structures/pnfine';
		# structurename = 'structures/pnrough';
		structurename = 'structures/pnhomo';
		# no PML required
		dxp = 1.0;
		dxn = 1.0;
		
		ndop  =  1e19;
		pdop  = -1e19;
		layers.append(['contact_0', 'GaAs', 0.0,2*dxp, dxp,  pdop]);
		layers.append(['left_p'   , 'GaAs', 0.0, 15.0, dxp,  pdop]);
		layers.append(['barrier1' , 'GaAs', 0.0,  7.0, dxp,  1e10]); 
		layers.append(['middle'   , 'GaAs', 0.0,  6.0, dxp,  1e10]); 
		layers.append(['barrier2' , 'GaAs', 0.0,  7.0, dxn,  1e10]);
		layers.append(['right_n'  , 'GaAs', 0.0, 15.0, dxn,  ndop]);
		layers.append(['contact_1', 'GaAs', 0.0,2*dxn, dxn,  ndop]);
	if IWCE:
		structurename = 'structures/iwce';
		if PML:
			structurename = structurename+'_pml';
		
		# syntax: name, material, xmole, length[nm], dx[nm], doping [cm-3]
		ndop  =  1e19;
		pdop  = -1e19;
		bmf   = 0.3;
		qwmf  = 0.0;
		bmat  = 'AlGaAs';
		qwmat = 'GaAs';
		if PML:
			# same discretization as NEGF grid
			layers.append(['PML_mx'   ,  bmat,  bmf, 15.0, 0.2,  pdop]);
			layers.append(['left_p'   ,  bmat,  bmf, 10.0, 0.5,  pdop]);
			layers.append(['left_i'   ,  bmat,  bmf,  5.0, 0.5,  1e10]);
			layers.append(['barrier1' ,  bmat,  bmf,  3.0, 0.5,  1e10]); 
			layers.append(['iface1a'  ,  bmat,  bmf,  0.2, 0.1,  1e10]); 
			layers.append(['iface1b'  , qwmat, qwmf,  0.2, 0.1,  1e10]); 
			layers.append(['mid_bound', qwmat, qwmf,  8.0, 0.5,  1e10]); 
			layers.append(['iface2b'  , qwmat, qwmf,  0.2, 0.1,  1e10]);
			layers.append(['iface2a'  ,  bmat,  bmf,  0.2, 0.1,  1e10]); 
			layers.append(['barrier2' ,  bmat,  bmf,  3.0, 1.0,  1e10]);
			layers.append(['right_i'  ,  bmat,  bmf,  5.0, 1.0,  1e10]);
			layers.append(['right_n'  ,  bmat,  bmf, 10.0, 1.0,  ndop]);
			layers.append(['PML_px'   ,  bmat,  bmf, 15.0, 0.2,  ndop]);
		else:
			layers.append(['contact_0',  bmat,  bmf,  1.0, 0.5,  pdop]);
			layers.append(['left_p'   ,  bmat,  bmf, 10.0, 0.5,  pdop]);
			layers.append(['left_i'   ,  bmat,  bmf,  5.0, 0.5,  1e10]);
			layers.append(['barrier1' ,  bmat,  bmf,  3.0, 0.5,  1e10]); 
			layers.append(['iface1a'  ,  bmat,  bmf,  0.2, 0.1,  1e10]); 
			layers.append(['iface1b'  , qwmat, qwmf,  0.2, 0.1,  1e10]); 
			layers.append(['mid_bound', qwmat, qwmf,  8.0, 0.5,  1e10]); 
			layers.append(['iface2b'  , qwmat, qwmf,  0.2, 0.1,  1e10]);
			layers.append(['iface2a'  ,  bmat,  bmf,  0.2, 0.1,  1e10]); 
			layers.append(['barrier2' ,  bmat,  bmf,  3.0, 1.0,  1e10]);
			layers.append(['right_i'  ,  bmat,  bmf,  5.0, 1.0,  1e10]);
			layers.append(['right_n'  ,  bmat,  bmf, 10.0, 1.0,  ndop]);
			layers.append(['contact_1',  bmat,  bmf,  2.0, 1.0,  ndop]);
	if NICHIA:
		structurename = 'structures/nichia';
		#structurename = 'structures/nichia_fine';
		if PML:
			structurename = structurename+'_pml';
		
		# syntax: name, material, xmole, length[nm], dx[nm], doping [cm-3]
		ndop  =  1e19;
		pdop  = -1e19;
		b1mf  = 0.0;
		qwmf  = 0.2;
		b2mf  = 0.2;
		b1mat = 'GaN';
		qwmat = 'InGaN';
		b2mat = 'AlGaN';
		dx_n = 0.4;
		dx_p = 0.25;
		if PML:
			# same discretization as NEGF grid
			layers.append(['PML_mx'   , b1mat, b1mf, 14.0, 0.2,  ndop]);
			layers.append(['left_n'   , b1mat, b1mf, 14.0, dx_n,  ndop]);
			#layers.append(['left_i'   , b1mat, b1mf,  5.0, dx_n,  1e10]);
			#layers.append(['barrier1' , b1mat, b1mf,  3.0, dx_n,  1e10]); 
			layers.append(['iface1a'  , b1mat, b1mf,  0.2, 0.05,  1e10]); 
			layers.append(['iface1b'  , qwmat, qwmf,  0.2, 0.05,  1e10]); 
			layers.append(['mid_bound', qwmat, qwmf,  3.0, dx_p,  1e10]); 
			layers.append(['iface2b'  , qwmat, qwmf,  0.2, 0.05,  1e10]);
			layers.append(['iface2a'  , b2mat, b2mf,  0.2, 0.05,  1e10]); 
			#layers.append(['barrier2' , b2mat, b2mf,  3.0, dx_p,  1e10]);
			#layers.append(['right_i'  , b2mat, b2mf,  5.0, dx_p,  1e10]);
			layers.append(['right_p'  , b2mat, b2mf, 10.0, dx_p,  pdop]);
			layers.append(['PML_px'   , b2mat, b2mf, 15.0, dx_p,  pdop]);
		else:
			layers.append(['contact_0', b1mat, b1mf,2*dx_n, dx_n,  ndop]);
			layers.append(['left_n'   , b1mat, b1mf,  14.0, dx_n,  ndop]);
			#layers.append(['left_i'   , b1mat, b1mf,   5.0, dx_n,  1e10]);
			#layers.append(['barrier1' , b1mat, b1mf,   3.0, dx_n,  1e10]); 
			layers.append(['iface1a'  , b1mat, b1mf,   0.2, 0.05,  1e10]); 
			layers.append(['iface1b'  , qwmat, qwmf,   0.2, 0.05,  1e10]); 
			layers.append(['mid_bound', qwmat, qwmf,   3.0, dx_p,  1e10]); 
			layers.append(['iface2b'  , qwmat, qwmf,   0.2, 0.05,  1e10]);
			layers.append(['iface2a'  , b2mat, b2mf,   0.2, 0.05,  1e10]); 
			#layers.append(['barrier2' , b2mat, b2mf,   3.0, dx_p,  1e10]);
			#layers.append(['right_i'  , b2mat, b2mf,   5.0, dx_p,  1e10]);
			layers.append(['right_p'  , b2mat, b2mf,  12.0, dx_p,  pdop]);
			layers.append(['contact_1', b2mat, b2mf,2*dx_p, dx_p,  pdop]);
	# =====================================================================
	# END OF USER INPUT
	# =====================================================================

	# ------------------------------------------------------------
	# create list of vertices, edges, elements, regions, contacts
	# also set up doping and xmole vector
	# ------------------------------------------------------------
	# initialize stuff
	conv = constants.convert_from_SI(units.length, 1e-9);
	x = -layers[0][3]; # leftmost x-coordinate!
	vertices = [];
	vcount = 0;
	edges = [];
	edcount = 0;
	elements = [];
	elcount = 0;
	regions = [];
	rcount = 0;
	contacts = [];
	ncontacts = 0;
	doping = DoubleVector();
	ndoping = DoubleVector();
	pdoping = DoubleVector();
	xmole = DoubleVector();
	
	cvar.logmsg.set_level(LOG_INFO);
	
	# do it!
	vertices.append(Vertex(vcount,x));
	xmole.push_back(layers[0][2]);
	dop = constants.convert_from_SI(units.density_3d, layers[0][5]*1e6);
	doping.push_back(dop);
	if dop>=0.0:
		ndoping.push_back(dop);
		pdoping.push_back(0.0);
	else:
		ndoping.push_back(0.0);
		pdoping.push_back(-dop);
		
	vcount += 1;
	for ll in range(len(layers)):
		r = Region(layers[ll][0]);
		r.set_index(rcount);
		r.set_material_name(layers[ll][1]);
		r.set_material_molefraction(layers[ll][2]);
		regions.append(r);
		rcount += 1;
		if ll==0:
			c = Contact('Source');
			c.set_index(0);
			contacts.append(c);
			ncontacts += 1;
		if ll==len(layers)-1:
			c = Contact('Drain');
			c.set_index(1);
			contacts.append(c);
			ncontacts += 1;
		
		width = layers[ll][3];
		dx = layers[ll][4];
		num_intervals = int(width / dx);
		cvar.logmsg.emit(LOG_INFO,"created region %s (molefraction %g) which will contain %d intervals." % (layers[ll][0],layers[ll][2],num_intervals));
		if abs(num_intervals - (width / dx)) > 1e-10:
			raise StandardError, 'layer '+str(ll)+': width ('+str(width)+') must be multiple of dx ('+str(dx)+'); width/dx='+str(width/dx);
		for ii in range(num_intervals):
			x += dx;
			vertices.append(Vertex(vcount,x * conv));
			vcount += 1;
			cvar.logmsg.emit(LOG_INFO,"created vertex %d (x=%.1fnm)" % (vcount,x));
			xmole.push_back(layers[ll][2]);
			dop = constants.convert_from_SI(units.density_3d, layers[ll][5]*1e6);
			doping.push_back(dop);
			if dop>=0.0:
				ndoping.push_back(dop);
				pdoping.push_back(0.0);
			else:
				ndoping.push_back(0.0);
				pdoping.push_back(-dop);
			
			edges.append(Edge(edcount, vertices[len(vertices)-1],vertices[len(vertices)-2]));
			edcount += 1;
			cvar.logmsg.emit(LOG_INFO_L2,'created edge %d' % edcount);
			el = Element(elcount, element_type.interval);
			el.add_vertex(vertices[len(vertices)-1]);
			el.add_vertex(vertices[len(vertices)-2]);
			el.add_edge(edges[len(edges)-1]);
			elements.append(el);
			elcount += 1;
			r.add_element(el);
			el.set_region(r);
			cvar.logmsg.emit(LOG_INFO_L2,"created element %d and added to region %s" % (elcount,r.get_name()));
			
			# add single vertices to the contacts
			# last vertex of first layer
			if ll==0 and ii==num_intervals-1:
				c.add_vertex(vertices[len(vertices)-1]);
			
			# first vertex of last layer
			if ll==len(layers)-1 and ii==0:
				c.add_vertex(vertices[len(vertices)-2]);
	
	# -------------------------------------------------------
	# create Geometry object, add everything and prepare
	# -------------------------------------------------------
	geom = Geometry(vcount, edcount, 0, elcount);
	for ii in range(vcount):
		geom.add_vertex(vertices[ii]);
	for ii in range(edcount):
		geom.add_edge(edges[ii]);
	for ii in range(elcount):
		geom.add_element(elements[ii]);
	for ii in range(rcount):
		geom.add_region(regions[ii]);
	for ii in range(ncontacts):
		geom.add_contact(contacts[ii]);
	#geom.set_num_dfise_elems(elcount);
	geom.set_num_dfise_elems(elcount + ncontacts);
	geom.set_num_dfise_regions(rcount + ncontacts);
	geom.prepare();
	geom.verify();
	
	# -------------------------------------------------------
	# write geometry to file
	# -------------------------------------------------------
	parser = InputParser();
	parser.write_dfise_grd(structurename, geom);
	
	# -------------------------------------------------------
	# set up doping equation and xmolefraction equation
	# and everything needed to write equations to file
	# -------------------------------------------------------
	doping_eqn = ArbitraryData(quantities.hole_density, doping);
	xmole_eqn = ArbitraryData(quantities.particles, xmole); # particles are unitless
	values = DoubleDoubleVector(); # const vector<const vector<double> *> 
	values.push_back(doping);
	values.push_back(ndoping);
	values.push_back(pdoping);
	values.push_back(xmole);
	datanames = StringVector();
	datanames.push_back('DopingConcentration');
	datanames.push_back('PhosphorusActiveConcentration');
	datanames.push_back('BoronActiveConcentration');
	datanames.push_back('xMoleFraction');
	locations = StringVector();
	locations.push_back('vertex');
	locations.push_back('vertex');
	locations.push_back('vertex');
	locations.push_back('vertex');
	unittypes = UnitVector();
	unittypes.push_back(units.density_3d);
	unittypes.push_back(units.density_3d);
	unittypes.push_back(units.density_3d);
	unittypes.push_back(units.unitless);
	
	values2 = parser.get_best_vector(values);
	parser.write_dfise_dat(structurename+'.dat', geom, values2, 4, datanames, locations, unittypes);
	
	cvar.logmsg.emit(LOG_INFO,"Thank you for using NEGF. Please report bugs and comments to steiger@iis.ee.ethz.ch.");
	cvar.logmsg.emit(LOG_INFO,"====================================================================================");
	cvar.mpi.terminate();
	
except RuntimeError, e:
	print "RuntimeError caught!"
	print e
	cvar.logmsg.emit_all(LOG_INFO, "Thread: %d" % cvar.mpi.get_rank());
	cvar.logmsg.emit_all(LOG_INFO, "Logfile: %s" % cvar.fnames.get_logfile());
	cvar.mpi.terminate();
	sys.exit();

except StandardError, e:
	print e
	cvar.logmsg.emit_all(LOG_INFO, "Thread: %d" % cvar.mpi.get_rank());
	cvar.logmsg.emit_all(LOG_INFO, "Logfile: %s" % cvar.fnames.get_logfile());
	cvar.mpi.terminate();
	sys.exit();
