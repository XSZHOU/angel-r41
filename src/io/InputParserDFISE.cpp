// included in InputParser.cpp except when NODFISE flag is set


// ---------------------------
// DF-ISE routines
// ---------------------------

/** read grid file filename and create a geometry object
 * @param grid file name
 * @return on success, geometry object containing nodes, elements etc. on failure, NULL will be returned
 * */
Geometry* InputParser::read_dfise_grd(const char* filename) const throw (Exception *)
try {
    string   file(filename);
    ifstream ftest;
    ostringstream ssout;
    file.append(".grd");

    ftest.open(file.c_str());
    if(ftest) {
        logmsg->emit(LOG_INFO,  "Reading geometry from file: %s", file.c_str());
        ftest.close();
    } else {
        logmsg->emit(LOG_ERROR,  "Error while trying to read %s",file.c_str());
        NEGF_EXCEPTION("cannot open gridfile");
        return 0;
    }

    DFISEGrdReader * df_geom = new DFISEGrdReader(file.c_str());

    // Checking input validity
    NEGF_ASSERT(df_geom->get_dimension()>=1 && df_geom->get_dimension()<=3, "invalid grid file with strange dimension");
    NEGF_ASSERT(df_geom->get_num_vertices() > 0, "invalid grid file with no vertices");
    NEGF_ASSERT(df_geom->get_num_elements() > 0, "invalid grid file with no elements");
    NEGF_ASSERT(df_geom->get_num_regions()  > 0, "invalid grid file with no regions");
    NEGF_ASSERT(!(df_geom->get_num_edges()==0 && df_geom->get_dimension()>1), "invalid grid file with no edges");
    NEGF_ASSERT(!(df_geom->get_dimension()==3 && df_geom->get_num_faces()==0), "invalid 3D-grid file with no faces");

    logmsg->emit(LOG_INFO,  "   DF-ISE geometry has %d vertices, %d elements and %d regions.",
            df_geom->get_num_vertices(), df_geom->get_num_elements(), df_geom->get_num_regions());

    // ------------------------------------------------------------------------------
    // determine how many 'real' elements there are
    // (contact elements have lower dimensionality and will be kicked out later on)
    // ------------------------------------------------------------------------------
    element_type::ElementType eltype;
    uint el_dim;
    uint valid_elems = 0;
    for(int eid = 0; eid < df_geom->get_num_elements(); eid++) {
        this->dfise_to_new( df_geom->get_element_type(eid), eltype, el_dim);
        if ((unsigned int)df_geom->get_dimension()==el_dim)
            valid_elems++;
    }
    logmsg->emit(LOG_INFO_L2,  "%d elements will be discarded (contact elems)",
                        df_geom->get_num_elements()-valid_elems, df_geom->get_num_elements());

    // -----------------------------------------------------
    // create the geometry object
    // -----------------------------------------------------
    Geometry * geometry = 0;
    if (df_geom->get_dimension()!=1) {
        geometry = new Geometry(  df_geom->get_num_vertices(), df_geom->get_num_edges(), df_geom->get_num_faces(), valid_elems);
    } else {    // special treatment of 1D because DF-ISE has no edges for it
        geometry = new Geometry(  df_geom->get_num_vertices(), valid_elems, df_geom->get_num_faces(), valid_elems);
    }
    geometry->set_dimension(df_geom->get_dimension());

    // -------------------------------------------------------------------
    // check for quasi-1D geometry NOT USED IN NEGF
    // -------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // add vertices
    // -------------------------------------------------------------------------

    logmsg->emit(LOG_INFO_L2, "Adding vertices...");

    for(int ii = 0; ii < df_geom->get_num_vertices(); ii++) {
        switch (geometry->get_dimension())
        {
        case 1:
            geometry->add_vertex( new Vertex(ii,
                constants::convert_from_SI(units::length, 1e-6*df_geom->get_vx(ii))) );
            break;
        case 2:
            geometry->add_vertex(new Vertex(ii,
                constants::convert_from_SI(units::length, 1e-6*df_geom->get_vx(ii)),
                constants::convert_from_SI(units::length, 1e-6*df_geom->get_vy(ii))) );
            break;
        case 3:
            geometry->add_vertex(new Vertex(ii,
                constants::convert_from_SI(units::length, 1e-6*df_geom->get_vx(ii)),
                constants::convert_from_SI(units::length, 1e-6*df_geom->get_vy(ii)),
                constants::convert_from_SI(units::length, 1e-6*df_geom->get_vz(ii))) );
        }
        geometry->get_vertex(ii)->set_index_external(ii);
    }

    // -------------------------------------------------------------------
    // add edges
    // -------------------------------------------------------------------
    int nvert = 0;
    Edge * edg;
    if (df_geom->get_dimension()==1) {
        NEGF_ASSERT(df_geom->get_num_edges()==0, "A 1D grid is not supposed to have edges.");
        // edges will be created later according to 1 element = 1 edge
    }
    else {
        logmsg->emit(LOG_INFO_L2, "Adding edges...");
    }
    for (int eid = 0; eid < df_geom->get_num_edges(); eid++)
    {
        const vector<uint> & vertices2 = df_geom->get_edge_vertices(eid);
        edg = new Edge(eid,geometry->get_vertex(vertices2[0]),geometry->get_vertex(vertices2[1]));
        edg->set_index_external(eid);
        geometry->add_edge(edg);
    }

    // -------------------------------------------------------------------
    // add faces and determine boundary vertices
    // -------------------------------------------------------------------
    nvert      = 0;
    int nedge  = 0;
    char location;
    int num_faces = df_geom->get_num_faces();
    int next  = 0;
    int count = 0;
    Face * fac;
    switch (df_geom->get_dimension())
    {
    case 1:
        NEGF_ASSERT(num_faces==0, "A 1D grid is not supposed to have faces");
        for (int elem_id = 0; elem_id < df_geom->get_num_elements(); elem_id++)
        {
            vector<uint> vertices; df_geom->get_element_vertices(elem_id, vertices);
            for (uint vid = 0; vid < vertices.size(); vid++)
            {
                uint idx = /*vid*/ geometry->get_vertex(vertices[vid])->get_index_global();
                location = df_geom->get_location(idx);  // 1D --> locations are defined for vertices
                if(location=='e') {
                    logmsg->emit(LOG_INFO_L2,"Vertex %d is external.", idx);
                    geometry->get_vertex(vertices[vid])->set_index_internal(-1);
                }
            }
        }
        break;
    case 2:
        NEGF_ASSERT(num_faces==0, "A 2D grid is not supposed to have faces");
        for (uint edge_id = 0; edge_id < (uint)df_geom->get_num_edges(); edge_id++)
        {
            location = df_geom->get_location(edge_id);  // 2D --> locations are defined for edges
            if(location=='e')
            {
                const vector<uint> & vertices = df_geom->get_edge_vertices(edge_id);
                for (uint jj = 0; jj < 2; jj++) {
                    geometry->get_vertex(vertices[jj])->set_index_internal(-1);
                }
            }
        }
        break;
    case 3:
        logmsg->init_progress_bar(LOG_INFO_L1,"Adding faces and determing boundary vertices...", num_faces);
        for(int fid = 0; fid < num_faces; fid++) {
            vector<uint> vertices; df_geom->get_face_vertices(fid, vertices);
            nvert = vertices.size();

            vector<uint> edges; df_geom->get_face_edges(fid, edges);
            nedge = edges.size();

            fac = new Face(fid, nvert, nedge);
            fac->set_index_external(fid);

            // add vertices to the face object and mark them if the face is exterior
            location = df_geom->get_location(fid);      // 3D --> locations are defined for faces
            // i = interior face in a region, f = face between regions, e = face on outside
            for(uint jj = 0; jj < (uint)nvert; jj++) {
                NEGF_ASSERT(vertices[jj] >= 0, "vertices[jj] >= 0");
                fac->add_vertex(geometry->get_vertex(vertices[jj]));
                if(location=='e')
                    (geometry->get_vertex(vertices[jj]))->set_index_internal(-1);
            }

            // add edges to the face object
            for(unsigned int jj = 0; jj < (unsigned int)nedge; jj++) {
                NEGF_ASSERT(edges[jj] >= 0, "edges[jj] >= 0");
                fac->add_edge(geometry->get_edge(edges[jj]));
            }

            geometry->add_face(fac);

            if(next == count++) next = logmsg->set_progress_bar(count, num_faces);
        }
        logmsg->end_progress_bar();
        break;
    default:
        NEGF_EXCEPTION("Strange dimensionality.");
        break;
    }

    // ----------------------------
    // add elements
    // ----------------------------
    nvert      = 0;
    nedge      = 0;
    int nface  = 0;
    int* new_indices = new int[ df_geom->get_num_elements() ]; // gives the new index of DF_ISE elemt ii (-1 if kicked out); needed for regions
    geometry->set_num_dfise_elems(df_geom->get_num_elements()); // only for mapping back results onto DF-ISE grid
    next       = 0;
    count      = 0;
    uint elem_idx = 0;
    logmsg->init_progress_bar(LOG_INFO_L1,"Adding elements...", df_geom->get_num_elements());
    for(int eid = 0; eid < df_geom->get_num_elements(); eid++)
    {
        // get element type and dimensionality
        int df_type = df_geom->get_element_type(eid);
        this->dfise_to_new( df_type, eltype, el_dim);

        vector<uint> vertices; df_geom->get_element_vertices(eid, vertices);
        nvert = vertices.size();

        vector<uint> edges; df_geom->get_element_edges(eid, edges);
        nedge = edges.size();

        //faces      = df_geom->get_ELallfaces(eid,nface); // creates new array on heap EXCEPT IF NFACE=0
        //NEGF_ASSERT(!(el_dim==3 && nface == 0), "nface = 0 for a 3D-element");
        vector<uint> faces; df_geom->get_element_faces(eid, faces);
        nface = faces.size();

        // handling of contact elements (in DF-ISE defined as elements of lower dimensionality)
        if ((unsigned int)df_geom->get_dimension() != el_dim)
        {
            if (el_dim != (unsigned int)df_geom->get_dimension()-1) {
                cout << "eid=" << eid << "eldim=" << el_dim << endl;
                NEGF_EXCEPTION("Don't know what to do with this element.");
            }

            // do not add the element to the geometry.
            new_indices[eid] = -1; // cout << "marked element "<<eid<<" as contact element. "<<endl;

            continue;
        }

        Element * elem;
        if (df_geom->get_dimension()==1)    // special treatment for 1D because DF-ISE interval elements don't have edges
        {
            NEGF_ASSERT(nedge==0 && nvert==2 && eltype==element_type::interval && nface==0, "1D elements must be intervals.");
            elem = new Element(elem_idx, eltype);
        } else {
            elem = new Element(elem_idx, eltype); // destructor is called by geom object
        }
        elem->set_index_external(eid);
        new_indices[eid] = elem_idx;

        // add vertices to the element object
        for(int jj = 0; jj < nvert; jj++) {
            NEGF_ASSERT(vertices[jj] >= 0, "vertices[jj] >= 0");
            elem->add_vertex(geometry->get_vertex((unsigned int)(vertices[jj])));
        }

        // add edges to the element object
        for(int jj = 0; jj < nedge; jj++) {
            if (edges[jj]<0)
                cout << "element "<<eid<<" edges["<<jj<<"]="<<edges[jj]<<endl;
            NEGF_ASSERT(edges[jj] >= 0, "edges[jj] >= 0");
            elem->add_edge(geometry->get_edge((unsigned int)(edges[jj])));
        }
        if (df_geom->get_dimension()==1)    // special treatment for 1D because DF-ISE interval elements don't have edges
        {
            edg = new Edge(elem_idx,elem->get_vertex(0),elem->get_vertex(1));
            edg->set_index_external(eid);
            elem->add_edge(edg);
            geometry->add_edge(edg);
        }

        // add faces to the element object
        for(int jj = 0; jj < nface; jj++) {
            NEGF_ASSERT(faces[jj] >= 0, "faces[jj] >= 0");
            elem->add_face(geometry->get_face((unsigned int)(faces[jj])));
        }

        geometry->add_element(elem);

        elem_idx++;
        if(next == count++)
            next = logmsg->set_progress_bar(count, df_geom->get_num_elements());
    }
    logmsg->end_progress_bar();
    NEGF_ASSERT(elem_idx==valid_elems, "Something went wrong during contact treatment.");

    // -------------------------------------------------------------------
    // add regions
    // -------------------------------------------------------------------
    logmsg->emit(LOG_INFO_L2, "Adding regions...");
    geometry->set_num_dfise_regions(df_geom->get_num_regions()); // only for mapping back results onto DF-ISE grid
    Region * region;
    for(int ii = 0; ii < df_geom->get_num_regions(); ii++) {
        // if the region's material is "Contact", don't add the region
        if (df_geom->get_region_material(ii)=="Contact") {
            continue;
        }

        // set region name and material name
        region = new Region(df_geom->get_region_name(ii));
        region->set_material_name(df_geom->get_region_material(ii).c_str());
        // assign region to elements
        for(uint jj = 0; jj < df_geom->get_num_region_elements(ii); jj++) {
            NEGF_ASSERT(df_geom->get_region_element(ii,jj) >= 0, "idx region element >= 0");
            int dfise_elem_idx = df_geom->get_region_element(ii,jj);
            if (new_indices[dfise_elem_idx] != -1) { // else the element was kicked out
                geometry->get_element(new_indices[dfise_elem_idx])->set_region(region);
            }
        }
        // assign elements to region
        for(uint jj = 0; jj < df_geom->get_num_region_elements(ii); jj++) {
            NEGF_ASSERT(df_geom->get_region_element(ii,jj) >= 0, "idx region element >= 0");
            uint dfise_elem_idx = df_geom->get_region_element(ii,jj);
            if (new_indices[dfise_elem_idx] != -1) { // else the element was kicked out
                region->add_element(geometry->get_element(new_indices[dfise_elem_idx]));
            }
        }
        // also add vertices to regions (necessary ONLY for reading DF-ISE .dat-files later on)
        vector<uint> region_verts; df_geom->get_region_vertices(ii, region_verts);
        nvert = region_verts.size();

        region->set_dfise_region_vertex_numbers(region_verts);
        geometry->add_region(region); // includes determination of region index and incrementing nregions
        logmsg->emit(LOG_INFO, "   added region '%s' (%d vertices)", region->get_name().c_str(), nvert);
    }
    delete [] new_indices; new_indices = 0;

    // -------------------------------------------------
    // add contacts (which are in DF-ISE also regions)
    // -----------------------------------------------
    logmsg->emit(LOG_INFO_L2, "Adding contacts...");
    Contact * contact = 0;
    for(int ii = 0; ii < df_geom->get_num_regions(); ii++) {
        if (df_geom->get_region_material(ii)!="Contact") {
            continue;
        }

        contact = new Contact(df_geom->get_region_name(ii));
        geometry->add_contact(contact); // numbering of contacts is done by add_contact

        // assign the contact to the vertices
        // IMPORTANT: we assume that the vertex index in our geometry index is the
        //            same as the vertex index in the DF-ISE file!!!
        vector<uint> vertices; df_geom->get_region_vertices(ii, vertices);
        nvert = vertices.size();

        for (int jj = 0; jj < nvert; jj++) {
            // geometry->get_vertex(vertices[jj])->set_contact(contact);
            contact->add_vertex(geometry->get_vertex(vertices[jj]));
            // add_vertex also assigns the contact and the index within the vertex list of the contact
            // to the contact
        }

        logmsg->emit(LOG_INFO_L2, "   added contact '%s' (%d vertices)",df_geom->get_region_name(ii).c_str(),nvert);
    }


    // ----------------------------------------------------------------------
    // add interfaces
    // NOT YET IMPLEMENTED!
    // ----------------------------------------------------------------------

    // ----------------------------------------------------------------------------------
    // check all elements. if their material corresponds to keyword "Gas", set it to zero
    // ----------------------------------------------------------------------------------
    const string gas("Gas");
    int count_gas = 0;
    int test      = 0;
    for(uint ii = 0; ii < geometry->get_num_elements(); ii++) {
        Element * elem = geometry->get_element(ii);
        if(elem->get_region()->get_material_name() ==   gas) {
            // set vertices to be dirichlet vertices
            for(uint jj = 0; jj < elem->get_num_vertices(); jj++) {
                (elem->get_vertex(jj))->set_index_internal(-1);
            }
            count_gas++;
        }
        test++;
    }

    if(count_gas > 0) {
        logmsg->emit(LOG_INFO_L1,  "%d elements assigned to material %s and their vertices will be ignored", count_gas, gas.c_str());
    }
    logmsg->emit(LOG_INFO_L2, "finished reading grid file");

    logmsg->emit(LOG_INFO_L2, "************* preparing geometry for calculation **************");
    geometry->prepare();

    geometry->verify();

    logmsg->emit(LOG_INFO_L2, "finished creating geometry");

    delete df_geom;
    return geometry;

} catch (       Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
  catch (dfise::Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }


void InputParser::get_fields(const char* filename, vector<string> & fieldnames,
                                vector<string> & locations) const throw (Exception *)
{//STACK_TRACE(
try {
    NEGF_ASSERT(filename!=0, "null pointer encountered.");
    string file(filename);
    file.append(".dat");

    DFISEDatReader dat_file(file.c_str());
    fieldnames = dat_file.get_dataset_names();          // const vector<string> & get_dataset_names() --> will copy array
    locations  = dat_file.get_locations();              // const vector<string> & get_locations() --> will copy array
    return;
} catch (       Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
  catch (dfise::Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
}
//);}


/** Check the existence of some field within a DF-ISE .dat-file */
bool InputParser::check_field_existence(const char* filename, const string & fieldname) const throw (Exception *)
{//STACK_TRACE(
try {
    NEGF_ASSERT(filename!=0 && fieldname!="", "wrong input to check_field_existence.");
    string file(filename);
    file.append(".dat");

    DFISEDatReader dat_file(file.c_str());
    if (dat_file.get_num_sets(fieldname) > 0) {
        return true;
    } else {
        return false;
    }
} catch (       Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
  catch (dfise::Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
}
//);}




/** Read a field from a DF-ISE .dat-file
 *  Since the geometry structure in this program is slightly different from DF-ISE,
 *  also the gridfile is needed to establish the correspondence.
 *  Furthermore it is assumed that a quantity is stored in separate fields for each region.
 *  @param filename the filename including path without .dat-ending
 *  @param gridname the corresponding DF-ISE grid including path without .grd-ending
 *  @param geom the Geometry object (used for correspondence of numbers)
 *  @param values the vector where the entries shall be stored (output)
 *  @param fieldname name of the field to read in
 *  @param loc the location of the field (vertex, element, ...)
 */
void InputParser::read_dfise_dat( const char* datfilename, const char* gridname,
                    const Geometry * const geom, vector<double> & values,
                    const string & fieldname, const string & loc) const throw (Exception *)
{
try {
    NEGF_ASSERT( geom!=0 && datfilename!=0 && fieldname!="" && loc!="", "wrong input to read_dfise_field.");
    values.clear();

    string file(datfilename);
    file.append(".dat");
    DFISEDatReader dat_file(file.c_str());

    // some checks
    NEGF_ASSERT(dat_file.get_dimension()==geom->get_dimension(),           "inconsistent dimension.");
    NEGF_ASSERT(dat_file.get_num_vertices()==geom->get_num_vertices(),     "inconsistent number of vertices.");
    if (geom->get_dimension()!=1) {
        NEGF_ASSERT(dat_file.get_num_edges()==geom->get_num_edges(),       "inconsistent number of edges.");
    }
    NEGF_ASSERT(dat_file.get_num_faces()==geom->get_num_faces(),           "inconsistent number of faces.");
    NEGF_ASSERT(dat_file.get_num_elements()==geom->get_num_dfise_elems(),  "inconsistent number of elements.");
    NEGF_ASSERT(dat_file.get_num_regions()==geom->get_num_dfise_regions(), "inconsistent number of regions.");

    // consistency checks
    uint num_sets = dat_file.get_num_sets(fieldname);
    NEGF_FASSERT(num_sets>0, "no dataset w/ name %s found in file %s", fieldname.c_str(), file.c_str());
    NEGF_ASSERT(dat_file.get_location(fieldname)=="vertex", "At the moment only fields defined on vertices can be read.");
    NEGF_ASSERT(dat_file.get_dataset_type(fieldname)=="scalar" && dat_file.get_dataset_dim(fieldname)==1, "At the moment only scalar fields can be read.");

    values.resize(geom->get_num_vertices(), 0.0);

    const vector< vector<double> > & dset_vals = dat_file.get_values(fieldname);
    NEGF_ASSERT(dset_vals.size()==num_sets, "inconsistency.");
    const vector< vector<string> > & dset_regions = dat_file.get_validities(fieldname);

    // read!
    if (num_sets==1) // assumption in this case: numbering is simply the global vertex numbering
    {
        for (uint ii = 0; ii < geom->get_num_vertices(); ii++)
            values[ii] = dset_vals[0][geom->get_vertex(ii)->get_index_external()];
    }
    else
    {
        // we assume that every set corresponds to a region
        // note that the field value at region interface vertices is random because the value
        // of the set with the highest set number will be taken.

        vector<double> dfise_numbering_values;
        dfise_numbering_values.resize(geom->get_num_vertices(), 0.0); // df_geom.NumVerts()==geom->get_num_vertices() was checked above
        for(uint ii = 0; ii < num_sets; ii++)
        {
            NEGF_ASSERT(dset_regions[ii].size()==1, "A set contained more than one region.");
            const vector<double> & region_vals = dset_vals[ii];

            // determine the DF-ISE region index of the set's region (-->ii2)
            // and the region index in NEGF (-->ii3)
            int ii3 = -1;
            for (uint jj = 0; jj < geom->get_num_regions(); jj++) {
                if (dset_regions[ii][0]==geom->get_region(jj)->get_name()) {
                    ii3 = jj;
                    break;
                }
            }
            NEGF_ASSERT( ii3 != -1, "No region with the same name could be found in NEGF.");

            NEGF_FASSERT(region_vals.size() == geom->get_region(ii3)->get_num_dfise_vertices(),
                    "Number of vertices in region \"%s\" is not consistent set->nb_values=%d, region->get_num_dfice_verts()=%d.",
                    geom->get_region(ii3)->get_name().c_str(), region_vals.size(), geom->get_region(ii3)->get_num_dfise_vertices());
            const vector<uint> & rverts = geom->get_region(ii3)->get_dfise_vertices();
            bool oxide_region = false;
            for (uint nn = 0; nn < Constants.oxide_names.size(); nn++) {
                if (geom->get_region(ii3)->get_material_name()==Constants.oxide_names[nn]) {    oxide_region = true; logmsg->emit(LOG_INFO,"REGION %s Is OXIDE",geom->get_region(ii3)->get_name().c_str()); break;  }
            }
            for(uint jj = 0; jj < region_vals.size(); jj++)
            {
                const vector<Region *> regs_nearby = geom->get_regions_near(geom->get_vertex(rverts[jj]));
                // handling of vertices at region interfaces
                // if the vertex is at an interface material-oxide/air, never take the oxide value
                bool skip_vertex = false;
                if (oxide_region && regs_nearby.size() > 1)
                {
                    bool non_oxide_region_around = false;
                    for (uint kk=0; kk<regs_nearby.size(); kk++) {
                        bool kk_oxide_region = false;
                        for (uint nn = 0; nn < Constants.oxide_names.size(); nn++) {
                            if (regs_nearby[kk]->get_material_name()==Constants.oxide_names[nn]) { kk_oxide_region = true; break; }
                        }
                        if (!kk_oxide_region) { non_oxide_region_around = true; break; }
                    }
                    if (non_oxide_region_around) { skip_vertex = true; }
                }
                if (!skip_vertex) {
                    dfise_numbering_values[rverts[jj]] = region_vals[jj];  // Assumption: Numbering here and in DF-ISE is the same
                } else {
                    logmsg->emit(LOG_INFO_L3,"Skipping value of field %s in region %s at vertex %d", fieldname.c_str(),
                            geom->get_region(ii3)->get_name().c_str(), rverts[jj]);
                }
            }
        }
        for (uint ii = 0; ii < geom->get_num_vertices(); ii++) {
            values[ii] = dfise_numbering_values[geom->get_vertex(ii)->get_index_external()];
        }
    }
} catch (       Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
  catch (dfise::Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
}


void InputParser::dfise_to_new(const int & dfise_type, element_type::ElementType & new_type, uint & el_dim) const
{STACK_TRACE(
    // conversion DFISE-element type --> new element type
    // DOUBLE-CHECK THIS CORRESPONDENCE!
    // for the correspondence compare the ENUM of Element.h with the ENUM of dfgeo.h
    switch (dfise_type)
    {
    case 0: // _df_point
        new_type = element_type::point;
        el_dim = 0;
        break;
    case 1: // _df_segment
        new_type = element_type::interval;
        el_dim = 1;
        break;
    case 2: // _df_triangle
        new_type = element_type::triangle;
        el_dim = 2;
        break;
    case 3: // _df_rectangle
        new_type = element_type::rectangle;
        el_dim = 2;
        break;
    case 4: // _df_polygon
        NEGF_EXCEPTION("An element type was read from file which is not supported.");
        break;
    case 5: // _df_tetrahedron
        new_type = element_type::tetrahedron;
        el_dim = 3;
        break;
    case 6: // _df_pyramid
        new_type = element_type::pyramid;
        el_dim = 3;
        break;
    case 7: // _df_prism
        new_type = element_type::prism;
        el_dim = 3;
        break;
    case 8: // _df_brick
        NEGF_EXCEPTION("An element type was read from file which is not supported.");
        break;
    case 9: // _df_tetrabrick
        new_type = element_type::tetrabrick;
        el_dim = 3;
        break;
    case 10: // _df_polyhedron
        NEGF_EXCEPTION("An element type was read from file which is not supported.");
        break;
    default:
        NEGF_EXCEPTION("An element type was read from file which is not supported.");
    }
    return;
);}


/** write (scalar) values to file in DF-ISE format
 *  @param filename     where to write
 *  @param geom         the Geometry object the quantity is defined on
 *  @param values_array a vector of pointers to the vectors containing the quantities to write.
 *  @param num_datasets how many datasets (fields) are written
 *  @param datanames    a vector of string containing the names of the datasets
 *  @param locations    a vector of strings containing where the quantities are located (vertex, element, ...)
 */
void InputParser::write_dfise_dat(const char* filename, const Geometry * geom,
                            const vector<const vector<double> *> values_array,
                            const uint num_datasets,
                            const vector<string> & datanames,
                            const vector<string> & locations,
                            const vector<units::UnitType>   & unittypes) const throw (Exception *)
{
    try {
    NEGF_ASSERT(filename!=0 && geom!=0, "wrong input to write_to_dfise_ile.");
    NEGF_ASSERT(datanames.size() == num_datasets, "size of name array must equal number of datasets.");
    NEGF_ASSERT(locations.size() == num_datasets, "size of location array must equal number of datasets.");
    NEGF_ASSERT(unittypes.size() == num_datasets, "size of units array must equal number of datasets.");
    NEGF_ASSERT(num_datasets > 0, "there is nothing to do!");

    DFISEDatWriter writer(geom->get_dimension(), geom->get_num_vertices(),
                        (geom->get_dimension()>1) ? geom->get_num_edges() : 0,
                        geom->get_num_faces(), geom->get_num_dfise_elems(),
                        geom->get_num_regions() + geom->get_num_contacts());

    vector<string> regionnames;
    for(uint ii=0; ii< geom->get_num_regions();  ii++) regionnames.push_back(geom->get_region(ii)->get_name());
    for(uint ii=0; ii< geom->get_num_contacts(); ii++) regionnames.push_back(geom->get_contact(ii)->get_name());
    for (uint dd = 0; dd < num_datasets; dd++)
    {
        const vector<double> & NEGF_values = *(values_array[dd]);

        vector<double> vals;
        vals.resize(NEGF_values.size(), 0.0);
        for (uint ii=0; ii<NEGF_values.size(); ii++) {
            vals[ii] = NEGF_values[ii] / constants::convert_from_SI(unittypes[dd], 1.0);
        }
        // check for not-a-number and and infinite entries
        for (uint ii = 0; ii < vals.size(); ii++) {
            if (isnan(vals[ii])) { logmsg->emit(LOG_WARN,  "Dataset %s: value %d is NaN. Visualizing the .dat-file may be a problem.",  datanames[dd].c_str(), ii); break; }
            if (isinf(vals[ii])) { logmsg->emit(LOG_WARN,  "Dataset %s: value %d is Inf. Visualizing the .dat-file may be a problem.",  datanames[dd].c_str(), ii); break; }
        }

        if (locations[dd]=="vertex") {
            NEGF_FASSERT(vals.size() == geom->get_num_vertices(),
                "size of value vector (%d) does not equal the number of vertices (%d).", vals.size(), geom->get_num_vertices());
        }
        else if (locations[dd]=="edge") { // no discarded edges
            NEGF_ASSERT(vals.size() == geom->get_num_edges(), "size of value vector does not equal the number of edges.");
            if (geom->get_dimension()==1) {
                logmsg->emit(LOG_WARN,  "Skipping output of field %s defined on 1D edges because DF-ISE does not support this.", datanames[dd].c_str());
                continue;
            }
        }
        else if (locations[dd]=="element") {
            NEGF_ASSERT( vals.size() == geom->get_num_elements(), "size of value vector does not equal the number of elements.");

            // hack the vector sol to re-include lost contact elements:
            vector<double> vals2 = vals;
            vals.assign(geom->get_num_dfise_elems(), 0.0);

            if (geom->get_dimension()==1) {
                //if (geom->get_num_dfise_elems()!=0) { logmsg->emit(LOG_WARN, "warning: we expect a DF-ISE 1D geometry not to have any elements."); }
                logmsg->emit(LOG_WARN, "Skipping output of field %s defined on 1D elements because DF-ISE does not support this.", datanames[dd].c_str());
            } else {
                for (uint ii = 0; ii < geom->get_num_elements(); ii++) {
                    uint idx = geom->get_element(ii)->get_index_external();
                    NEGF_ASSERT(idx < geom->get_num_dfise_elems(), "invalid external index.");
                    NEGF_ASSERT(vals[idx]==0.0, "An external index occurred twice.");
                    vals[idx] = vals2[ii];
                }
            }
        }
        else { NEGF_EXCEPTION("Unsupported location of quantity on grid."); }

        writer.add_dataset(datanames[dd], datanames[dd], "scalar", 1, locations[dd], regionnames, vals);
    }
    writer.write(filename);
} catch (       Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
  catch (dfise::Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
}


void InputParser::read_dfise_plt(const char* filename,
                                vector<string> & fieldnames,
                                vector< vector<double> > & values) throw (Exception *)
{
try {
    DFISEPltReader reader(filename);
    fieldnames = reader.get_dataset_names();
    values.resize(fieldnames.size());
    for (uint ii=0; ii<fieldnames.size(); ii++) {
        values[ii] = reader.get_data(fieldnames[ii]);
    }
} catch (       Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
  catch (dfise::Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
}


/** Write DF-ISE .plt-file (see DFISE/include/dfplt.h)
 *  the data of dataset i, value j must be stored in values[j*num_datasets+i]. */
void InputParser::write_dfise_plt(  const char* filename,
                                    const uint num_values,
                                    const uint num_datasets,
                                    const vector<string> & datanames,
                                    const vector<units::UnitType> & unittypes,
                                    const vector<double> & values) const throw (Exception *)
{
try {
    NEGF_ASSERT(num_values>0 && num_datasets>0, "there is nothing to be done.");
    NEGF_ASSERT(filename!=NULL, "encountered null  pointer.");
    NEGF_ASSERT(datanames.size()==num_datasets, "num_datasets does not coincide with number of names.");
    NEGF_ASSERT(values.size()==num_values*num_datasets, "values.size()!=num_values*num_datasets");
    NEGF_ASSERT(unittypes.size()==datanames.size(), "sizes of unittype and dataname vectors do not coincide.");

    vector<double> tmp_vals; tmp_vals.resize(num_values);
    vector< vector<double> > vals; vals.resize(num_datasets, tmp_vals);
    for (uint ii=0; ii<num_datasets; ii++) {
        for (uint jj=0; jj<num_values; jj++) {
            vals[ii][jj] = values[jj*num_datasets+ii] / constants::convert_from_SI(unittypes[ii], 1.0);
        }
    }

    DFISEPltWriter writer;
    for (uint ii=0; ii<num_datasets; ii++) {
        writer.add_data(datanames[ii], vals[ii]);
    }
    writer.write(filename);
} catch (       Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
  catch (dfise::Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
}


void InputParser::write_dfise_grd(const char * filename, const Geometry * const grid) const throw (Exception *)
{//STACK_TRACE(
try {
    DFISEGrdWriter mygeo(grid->get_dimension());

    // ---------------------------------------------------
    // vertices
    // --------------------------------------------------
    logmsg->emit(LOG_INFO_L1, "Preparing vertices");
    mygeo.set_num_vertices(grid->get_num_vertices());
    const double um = constants::convert_from_SI(units::length, 1e-6); // 1 um in NEGF units
    vector<double> coord; coord.resize(grid->get_num_vertices());
    for (uint ii = 0; ii < grid->get_num_vertices(); ii++) {
        coord[ii]=grid->get_vertex(ii)->get_coordinate(0) / um;
    }
    mygeo.set_x_coords(coord);
    if (grid->get_dimension() >= 2) {
        for (uint ii = 0; ii < grid->get_num_vertices(); ii++) {
            coord[ii]=grid->get_vertex(ii)->get_coordinate(1) / um;
        }
        mygeo.set_y_coords(coord);
    }
    if (grid->get_dimension() >= 3) {
        for (uint ii = 0; ii < grid->get_num_vertices(); ii++) {
            coord[ii]=grid->get_vertex(ii)->get_coordinate(2) / um;
        }
        mygeo.set_z_coords(coord);
    }

    // ----------------------------------------------------
    // edges
    // ----------------------------------------------------
    logmsg->emit(LOG_INFO_L1, "Preparing edges");
    if (grid->get_dimension() > 1) {
        mygeo.set_num_edges(grid->get_num_edges());
        for(uint ii=0; ii<grid->get_num_edges(); ii++) {
            mygeo.set_edge_vertices(ii, grid->get_edge(ii)->get_lower_vertex()->get_index_global(), grid->get_edge(ii)->get_upper_vertex()->get_index_global());
        }
    } else {
        mygeo.set_num_edges(0);
    }

    // ----------------------------------------------------
    // faces
    // ----------------------------------------------------
    logmsg->emit(LOG_INFO_L1, "Preparing faces");
    mygeo.set_num_faces(grid->get_num_faces());
    for(uint ii=0; ii<grid->get_num_faces(); ii++) {
        vector<int> face_edges;
        Face * face = grid->get_face(ii);

        Edge * e0; Edge * e1; Edge * e2; Edge * e3;
        switch (face->get_num_edges())
        {
        case 3: // triangle
            // find right order of edges
            // order: e0=(v0,v1); e1=(v1,v2); e2=(v2,v0);   vi = vertex i of the element
            e0 = this->get_edge_containing_specific_vertices(face, 0, 1);
            e1 = this->get_edge_containing_specific_vertices(face, 1, 2);
            e2 = this->get_edge_containing_specific_vertices(face, 2, 0);

            // write indices (possibly negative!)
            face_edges.resize(3);
            face_edges[0] = this->get_dfise_ordered_index(e0, face->get_vertex(0));
            face_edges[1] = this->get_dfise_ordered_index(e1, face->get_vertex(1));
            face_edges[2] = this->get_dfise_ordered_index(e2, face->get_vertex(2));
            break;
        case 4: // rectangle
            // find right order of edges
            // order: e0=(v0,v1); e1=(v1,v2); e2=(v2,v0);   vi = vertex i of the element
            e0 = this->get_edge_containing_specific_vertices(face, 0, 1);
            e1 = this->get_edge_containing_specific_vertices(face, 1, 2);
            e2 = this->get_edge_containing_specific_vertices(face, 2, 3);
            e3 = this->get_edge_containing_specific_vertices(face, 3, 0);

            // write indices (possibly negative!)
            face_edges.resize(4);
            face_edges[0] = this->get_dfise_ordered_index(e0, face->get_vertex(0));
            face_edges[1] = this->get_dfise_ordered_index(e1, face->get_vertex(1));
            face_edges[2] = this->get_dfise_ordered_index(e2, face->get_vertex(2));
            face_edges[3] = this->get_dfise_ordered_index(e3, face->get_vertex(3));
            break;
        default:
            NEGF_EXCEPTION("Faces other than triangles/rectangles are not supported.");
        }

        mygeo.set_face_edges(ii,face_edges);
    }


    // ----------------------------------------------------
    // elements
    // ----------------------------------------------------
    logmsg->emit(LOG_INFO_L1, "Preparing elements");

    // determine number of contact elements
    uint num_contact_elems = 0;
    switch (grid->get_dimension())
    {
    case 1:
        // search points if they are at a contact
        // NEGF special: exclude points which are not connected to non-contact points
        for (uint ii = 0; ii < grid->get_num_vertices(); ii++) {
            if (grid->get_vertex(ii)->is_at_contact()) {
                const vector<Edge *> near_edges = grid->get_edges_near(grid->get_vertex(ii));
                bool is_interface = false;
                for (uint jj=0; jj<near_edges.size(); jj++) {
                    Vertex * v = (near_edges[jj]->get_lower_vertex()==grid->get_vertex(ii)) ? near_edges[jj]->get_upper_vertex() : near_edges[jj]->get_lower_vertex();
                    if (!v->is_at_contact()) {
                        is_interface = true;
                        break;
                    }
                }
                if (is_interface) {
                    num_contact_elems++;
                }
            }
        }
        break;
    case 2:
        // search edges if they are at contact
        for (uint ii = 0; ii < grid->get_num_edges(); ii++) {
            if (   grid->get_edge(ii)->get_lower_vertex()->is_at_contact()
                && grid->get_edge(ii)->get_upper_vertex()->is_at_contact())  {
                num_contact_elems++;
            }
        }
        break;
    case 3:
        // search faces if they are at contact
        for (uint ii = 0; ii < grid->get_num_faces(); ii++) {
            bool contact_face = true;
            for (uint jj = 0; jj < grid->get_face(ii)->get_num_vertices(); jj++) {
                if (!grid->get_face(ii)->get_vertex(jj)->is_at_contact()) {
                    contact_face = false;
                    break;
                }
            }
            if (contact_face) {
                num_contact_elems++;
            }
        }
        break;
    default: NEGF_EXCEPTION("Strange grid dimension. You should not have got here.");
    }

    // determine total number of elements
    mygeo.set_num_elements(grid->get_num_elements() + num_contact_elems);

    // fill ordinary elements
    Edge * e0; Edge * e1; Edge * e2; Edge * e3;
    for(uint ii=0; ii < grid->get_num_elements(); ii++) {
        vector<int> elem_entities;

        Element * elem = grid->get_element(ii);
        switch (grid->get_dimension())
        {
        case 1:
            NEGF_ASSERT(elem->get_type()==element_type::interval, "wrong 1D element.");
            elem_entities.resize(3);

            elem_entities[0] = 1; // interval
            elem_entities[1] = elem->get_vertex(0)->get_index_global();
            elem_entities[2] = elem->get_vertex(1)->get_index_global();
            break;
        case 2:
            elem_entities.resize(elem->get_num_edges() + 1);
            switch(elem->get_type()) {
                case element_type::triangle:
                    elem_entities[0] = 2;

                    // find right order of edges
                    // order: e0=(v0,v1); e1=(v1,v2); e2=(v2,v0);   vi = vertex i of the element
                    e0 = this->get_edge_containing_specific_vertices(elem, 0, 1);
                    e1 = this->get_edge_containing_specific_vertices(elem, 1, 2);
                    e2 = this->get_edge_containing_specific_vertices(elem, 2, 0);

                    // write indices (possibly negative!)
                    elem_entities[1] = this->get_dfise_ordered_index(e0, elem->get_vertex(0));
                    elem_entities[2] = this->get_dfise_ordered_index(e1, elem->get_vertex(1));
                    elem_entities[3] = this->get_dfise_ordered_index(e2, elem->get_vertex(2));

                    break;
                case element_type::rectangle:
                    elem_entities[0] = 3;

                    // find right order of edges
                    // order: e0=(v0,v1); e1=(v1,v2); e2=(v2,v3);  e3=(v3,v0);  vi = vertex i of the element
                    e0 = this->get_edge_containing_specific_vertices(elem, 0, 1);
                    e1 = this->get_edge_containing_specific_vertices(elem, 1, 2);
                    e2 = this->get_edge_containing_specific_vertices(elem, 2, 3);
                    e3 = this->get_edge_containing_specific_vertices(elem, 3, 0);

                    // write indices (possibly negative!)
                    elem_entities[1] = this->get_dfise_ordered_index(e0, elem->get_vertex(0));
                    elem_entities[2] = this->get_dfise_ordered_index(e1, elem->get_vertex(1));
                    elem_entities[3] = this->get_dfise_ordered_index(e2, elem->get_vertex(2));
                    elem_entities[4] = this->get_dfise_ordered_index(e3, elem->get_vertex(3));

                    break;
                default:
                    NEGF_EXCEPTION("wrong 2D element.");
                    break;
            }
            break;
        // in 3D the whole ordering thing is not done (left to the reader as exercise)
        case 3:
            elem_entities.resize(elem->get_num_faces() + 1);
            switch(elem->get_type()) {
                case element_type::tetrahedron:
                    elem_entities[0] = 5;
                    break;
                case element_type::pyramid:
                    elem_entities[0] = 6;
                    break;
                case element_type::prism:
                    elem_entities[0] = 7;
                    break;
                case element_type::tetrabrick:
                    elem_entities[0] = 9;
                    break;
                default:
                    NEGF_EXCEPTION("wrong 3D element.");
                    break;
            }
            for(uint jj=0; jj<elem->get_num_faces(); jj++) {
                elem_entities[jj+1] = elem->get_face(jj)->get_index_global();
            }
            break;
        default: NEGF_EXCEPTION("Wrong dimensionality.");
        }

        mygeo.set_element_entities(ii, elem_entities);
    }
    // fill contact elements
    uint count = 0;
    vector< vector<uint> > contact_elem_numbers;
    contact_elem_numbers.resize(grid->get_num_contacts());
    vector<int> elem_entities;
    switch (grid->get_dimension())
    {
    case 1:
        // go through contact vertices
        for (uint ii = 0; ii < grid->get_num_vertices(); ii++) {
            if (grid->get_vertex(ii)->is_at_contact()) {
                const vector<Edge *> near_edges = grid->get_edges_near(grid->get_vertex(ii));
                bool is_interface = false;
                for (uint jj=0; jj<near_edges.size(); jj++) {
                    Vertex * v = (near_edges[jj]->get_lower_vertex()==grid->get_vertex(ii)) ? near_edges[jj]->get_upper_vertex() : near_edges[jj]->get_lower_vertex();
                    if (!v->is_at_contact()) {
                        is_interface = true;
                        break;
                    }
                }
                if (is_interface) {
                    elem_entities.resize(2);
                    elem_entities[0] = 0;   // point
                    elem_entities[1] = grid->get_vertex(ii)->get_index_global();
                    mygeo.set_element_entities(grid->get_num_elements()+count, elem_entities);

                    contact_elem_numbers[grid->get_vertex(ii)->get_contact()->get_index()].push_back(grid->get_num_elements()+count);
                    count++;
                }
            }
        }
        break;
    case 2:
        // go through contact edges
        for (uint ii = 0; ii < grid->get_num_edges(); ii++) {
            if (   grid->get_edge(ii)->get_lower_vertex()->is_at_contact()
                && grid->get_edge(ii)->get_upper_vertex()->is_at_contact())
            {
                elem_entities.resize(3);
                elem_entities[0] = 1;   // interval
                elem_entities[1] = grid->get_edge(ii)->get_lower_vertex()->get_index_global();
                elem_entities[2] = grid->get_edge(ii)->get_upper_vertex()->get_index_global();
                mygeo.set_element_entities(grid->get_num_elements()+count, elem_entities);

                contact_elem_numbers[grid->get_edge(ii)->get_lower_vertex()->get_contact()->get_index()]
                    .push_back(grid->get_num_elements()+count);
                count++;
            }
        }
        break;
    case 3:
        // go through contact faces
        for (uint ii = 0; ii < grid->get_num_faces(); ii++) {
            Face * face = grid->get_face(ii);

            bool contact_face = true;
            for (uint jj = 0; jj < face->get_num_vertices(); jj++) {
                if (!face->get_vertex(jj)->is_at_contact()) {
                    contact_face = false;
                    break;
                }
            }
            if (contact_face)
            {
                elem_entities.resize(face->get_num_edges() + 1);
                switch(face->get_num_vertices()) {
                case 3:
                    elem_entities[0] = 2;   // triangle

                    // find right order of edges
                    // order: e0=(v0,v1); e1=(v1,v2); e2=(v2,v0);   vi = vertex i of the element
                    e0 = this->get_edge_containing_specific_vertices(face, 0, 1);
                    e1 = this->get_edge_containing_specific_vertices(face, 1, 2);
                    e2 = this->get_edge_containing_specific_vertices(face, 2, 0);

                    // write indices (possibly negative!)
                    elem_entities[1] = this->get_dfise_ordered_index(e0, face->get_vertex(0));
                    elem_entities[2] = this->get_dfise_ordered_index(e1, face->get_vertex(1));
                    elem_entities[3] = this->get_dfise_ordered_index(e2, face->get_vertex(2));

                    break;
                case 4:
                    elem_entities[0] = 3;   // rectangle

                    // find right order of edges
                    // order: e0=(v0,v1); e1=(v1,v2); e2=(v2,v0);   vi = vertex i of the element
                    e0 = this->get_edge_containing_specific_vertices(face, 0, 1);
                    e1 = this->get_edge_containing_specific_vertices(face, 1, 2);
                    e2 = this->get_edge_containing_specific_vertices(face, 2, 3);
                    e3 = this->get_edge_containing_specific_vertices(face, 3, 0);

                    // write indices (possibly negative!)
                    elem_entities[1] = this->get_dfise_ordered_index(e0, face->get_vertex(0));
                    elem_entities[2] = this->get_dfise_ordered_index(e1, face->get_vertex(1));
                    elem_entities[3] = this->get_dfise_ordered_index(e2, face->get_vertex(2));
                    elem_entities[4] = this->get_dfise_ordered_index(e2, face->get_vertex(3));
                    break;
                default:
                    NEGF_EXCEPTION("Contact faces other than triangles and rectangles are not implemented.");
                    break;
                }
                mygeo.set_element_entities(grid->get_num_elements()+count, elem_entities);

                contact_elem_numbers[face->get_vertex(0)->get_contact()->get_index()].push_back(grid->get_num_elements()+count);
                count++;
            }
        }
        break;
    default: NEGF_EXCEPTION("You should not have got here.");
    }
    NEGF_ASSERT(num_contact_elems == count, "inconsistency.");


    // ----------------------------------------------------
    // regions
    // ----------------------------------------------------
    logmsg->emit(LOG_INFO_L1, "Preparing contacts");

    mygeo.set_num_regions(grid->get_num_regions() + grid->get_num_contacts());

    // determine region and material names
    logmsg->emit(LOG_INFO_L1, "   region and material names");
    for (uint ii = 0; ii < grid->get_num_regions(); ii++) {
        mygeo.set_region_name    (ii, grid->get_region(ii)->get_name());
        mygeo.set_region_material(ii, grid->get_region(ii)->get_material_name());
    }
    for (uint ii = 0; ii < grid->get_num_contacts(); ii++) {
        mygeo.set_region_name    (grid->get_num_regions()+ii, grid->get_contact(ii)->get_name());
        mygeo.set_region_material(grid->get_num_regions()+ii, "Contact");
    }

    // determine region elements
    logmsg->emit(LOG_INFO_L1, "   region elements");
    vector< vector<uint> > elem_numbers;
    elem_numbers.resize(grid->get_num_regions());
    for (uint ii = 0; ii < grid->get_num_elements(); ii++) {
        elem_numbers[grid->get_element(ii)->get_region()->get_index()].push_back(grid->get_element(ii)->get_index_global());
    }
    for (uint ii = 0; ii < grid->get_num_regions(); ii++) {
        mygeo.set_region_elements(ii, elem_numbers[ii]);
    }

    // determine contact elements
    logmsg->emit(LOG_INFO_L1, "   contact elements");
    for (uint ii = 0; ii < grid->get_num_contacts(); ii++) {
        mygeo.set_region_elements(grid->get_num_regions()+ii, contact_elem_numbers[ii]);
    }

    // -----------------------------------------------------
    // locations
    // -----------------------------------------------------
    logmsg->emit(LOG_INFO_L1, "Preparing locations");
    vector<char> locations;
    switch(grid->get_dimension())
    {
    case 1:
        locations.resize(grid->get_num_vertices());
        for (uint ii = 0; ii < grid->get_num_vertices(); ii++)
        {
            if (grid->get_vertex(ii)->is_at_contact())
                locations[ii] = 'e';
            else if (grid->get_regions_near(grid->get_vertex(ii)).size() > 1)
                locations[ii] = 'f';
            else
                locations[ii] = 'i';
        }
        break;
    case 2:
        locations.resize(grid->get_num_edges());
        for (uint ii = 0; ii < grid->get_num_edges(); ii++)
        {
            if ( grid->get_elems_near(grid->get_edge(ii)).size() <=1 ) {
                locations[ii] = 'e';
            } else {
                bool interface = false;
                const vector<Element *> elems_near_edge = grid->get_elems_near(grid->get_edge(ii));
                Region * reg = elems_near_edge[0]->get_region();
                for (uint jj = 1; jj < elems_near_edge.size(); jj++) {
                    if (elems_near_edge[jj]->get_region()!=reg) {
                        interface = true;
                        break;
                    }
                }
                if (interface) {
                    locations[ii] = 'f';
                } else {
                    locations[ii] = 'i';
                }
            }
        }
        break;
    case 3:
        locations.resize(grid->get_num_faces());
        NEGF_ASSERT(grid->is_prepared(), "prepare grid first.");
        for (uint ii = 0; ii < grid->get_num_faces(); ii++)
        {
            Face * face = grid->get_face(ii);

            Element * first_elem    = 0;
            Region *  first_region  = 0;
            Element * second_elem   = 0;
            Region *  second_region = 0;
            for (uint jj = 0; jj < face->get_num_edges(); jj++)
            {
                // look if the adjacent elements have the face in their list
                const vector<Element *> & elems_near_edge = grid->get_elems_near(face->get_edge(jj));
                for (uint kk = 0; kk < elems_near_edge.size(); kk++)
                {
                    for (uint ll = 0; ll < elems_near_edge[kk]->get_num_faces(); ll++)
                    {
                        if (elems_near_edge[kk]->get_face(ll)==face)
                        {
                            if (first_elem == 0) {
                                first_elem   = elems_near_edge[kk];
                                first_region = elems_near_edge[kk]->get_region();
                                second_elem   = first_elem;
                                second_region = first_region;
                            } else if (first_elem != elems_near_edge[kk]) {
                                second_elem = elems_near_edge[kk];
                                if (second_elem->get_region() != first_region) {
                                    second_region = second_elem->get_region();
                                    break;  // now we're finished
                                }
                            }
                        }
                    }
                    if (first_region != second_region)  break;
                }
                if (first_region != second_region)  break;
            }
            if (first_elem==0) {
                for (uint jj = 0; jj < face->get_num_vertices(); jj++) {
                    Vertex * v = face->get_vertex(jj);
                    logmsg->emit(LOG_ERROR,"   face %d vertex %d is vertex %d(%e,%e,%e)",
                            ii,jj,v->get_index_global(),v->get_coordinate(0),v->get_coordinate(1),v->get_coordinate(2));
                }
                for (uint jj = 0; jj < face->get_num_edges(); jj++) {
                    logmsg->emit(LOG_ERROR,"   face %d edge %d has %d adjacent elems.",
                        ii,jj, grid->get_elems_near(face->get_edge(jj)).size());
                }
                NEGF_FEXCEPTION("no element adjacent to face %d was found.",ii);
            }

            if (second_elem==first_elem)
                locations[ii] = 'e';
            if (second_elem!=first_elem && second_region==first_region)
                locations[ii] = 'i';
            if (second_elem!=first_elem && second_region!=first_region)
                locations[ii] = 'f';
        }
        break;
    default: NEGF_EXCEPTION("Strange dimension."); break;
    }
    mygeo.set_locations(locations);

    logmsg->emit(LOG_INFO_L1, "Writing to file");
    string grdfile(filename);
    grdfile.append(".grd");
    mygeo.write(grdfile.c_str());

    logmsg->emit(LOG_INFO_L1, "Done!");
} catch (       Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
  catch (dfise::Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
}


uint InputParser::get_dfise_ordered_index(const Edge * edge, const Vertex * first_vertex) const
{STACK_TRACE(
    if (edge->get_lower_vertex()==first_vertex) {
        return edge->get_index_global();
    } else {
        NEGF_ASSERT(edge->get_upper_vertex()==first_vertex, "given vertex is not connected to edge.");
        return -edge->get_index_global() - 1;
    }
);}




/** prepare molefractions w/ .dat-file and existing geometry */
void InputParser::prepare_molefractions(const Geometry * const geom, const string & filename) const
{STACK_TRACE(
    NEGF_ASSERT(geom!=NULL, "null pointer encountered.");

    logmsg->emit(LOG_INFO_L2,"    Checking if ternary materials are present...");
    // check if a ternary material is present at all
    bool ternary_material_found = this->check_if_ternary_materials_exist(geom);
    if (!ternary_material_found) // there is nothing to do
        return;

    // read in field from file
    logmsg->emit(LOG_INFO_L2,"    Trying to read field \"xMoleFraction\" from %s...",filename.c_str());
    bool field_found = this->check_field_existence(filename.c_str(), "xMoleFraction");
    NEGF_ASSERT(field_found, "Field \"xMoleFraction\" not found even though a ternary material was encountered.");

    vector<double> tempvec;
    this->read_dfise_dat(filename.c_str(), filename.c_str(), geom, tempvec, "xMoleFraction", "vertex");
    NEGF_ASSERT(tempvec.size()==geom->get_num_vertices(), "error while reading xMoleFraction.");

    // get for each region the value of region-interior points
    logmsg->emit(LOG_INFO_L2,"    Reducing vertex-based quantities to region-based quantities...");
    vector<double> x_values = this->get_regionwise_constant_xmole(geom, tempvec);

    // assign molefractions to each region
    logmsg->emit(LOG_INFO_L2,"    Assigning mole fractions...");
    this->assign_xmole_to_regions(geom, x_values);
);}



