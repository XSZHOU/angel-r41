// included in OutputData.cpp


void OutputData::write_dat() const throw (Exception *)
{STACK_TRACE(
     if(this->dat_equations.size()==0) {
        // There is nothing to write.
        return;
     }
    NEGF_ASSERT(this->dat_equations.size()==this->dat_unittypes.size(), "something went wrong.");

    // prepare things needed by DF-ISE output routines
    uint                            num_datasets = 0;
    vector<const vector<double> *>  value_array;
    vector<string>                  datanames;
    vector<string>                  locations;
    vector<units::UnitType>         types;
    vector< vector<double> * >      created_fields;
    for (uint ii = 0; ii < this->dat_equations.size(); ii++)
    {
        Equation * eqn = dat_equations[ii];

        // some checks
        NEGF_ASSERT(eqn!=0, "Found null pointer.");
        NEGF_ASSERT(eqn->get_name().size()>0, "Must hand over a name.");
        if ((eqn->get_values()).size()==0) {
            eqn->compute_values(eqn->get_timestamp());
        }
        NEGF_FASSERT((eqn->get_values()).size()==eqn->get_num_variables(),
                    "Number of variable values of eqn \"%s\" does not coincide with length of value vector.",
                    eqn->get_name().c_str());

        num_datasets++;
        value_array.push_back(&(eqn->get_values()));
        datanames.push_back(eqn->get_name());
        types.push_back(dat_unittypes[ii]);

        // --------------------------------------------------------------------------------------------
        // find out location
        // Note: in 1D #edges=#elems ALWAYS!
        // DF-ISE usually does not have edges in 1D, we created them during the geometry setup process
        // therefore it is better to store the variables on the elements
        // --------------------------------------------------------------------------------------------
        uint D = grid->get_dimension();
        NEGF_FASSERT( grid->get_num_vertices()!=grid->get_num_edges(),
                        "cannot determine location of a quantity because #verts=#edges=%d", grid->get_num_edges());
        NEGF_FASSERT( grid->get_num_vertices()!=grid->get_num_elements(),
                        "cannot determine location of a quantity because #verts=#elems=%d", grid->get_num_vertices());
        if (D==1) {
            if (eqn->get_num_variables()==this->grid->get_num_vertices())
                locations.push_back("vertex");
            else if (eqn->get_num_variables()==this->grid->get_num_elements())
                locations.push_back("element");
            else
                NEGF_EXCEPTION("Equation must be defined on vertices, edges or elements. Other possibility: Eqn is for a different grid.");
        } else {
            NEGF_FASSERT( grid->get_num_edges()   !=grid->get_num_elements(),
                        "cannot determine location of a quantity because #edges=#elems=%d", grid->get_num_edges());
            if (eqn->get_num_variables()==this->grid->get_num_vertices()) {
                locations.push_back("vertex");
            } else if (eqn->get_num_variables()==this->grid->get_num_edges()) {
                // locations.push_back("edge");
                NEGF_EXCEPTION("Although DF-ISE supports location=edge, TECPLOT is not able to display them.");
            } else if (eqn->get_num_variables()==this->grid->get_num_elements()) {
                locations.push_back("element");
            } else if (eqn->get_num_variables()==this->grid->get_num_vertices()*D) {
                // special treatment of vector fields
                // assume component d of vertex i to be stored in Di+d, D=grid-dimension
                locations.push_back("vertex");

                // redefine current field to be the modulus
                vector<double> * field_modulus = new vector<double>;
                created_fields.push_back(field_modulus);
                field_modulus->resize(grid->get_num_vertices(), 0.0);
                for (uint jj = 0; jj < grid->get_num_vertices(); jj++) {
                    for (uint dd = 0; dd < D; dd++) {
                        (*field_modulus)[jj] += eqn->get_value(D*jj+dd) * eqn->get_value(D*jj+dd);
                    }
                    (*field_modulus)[jj] = std::sqrt((*field_modulus)[jj]);
                }
                value_array[value_array.size()-1] = field_modulus;
                datanames[datanames.size()-1] = eqn->get_name() + "_mod";

                // add the separate components as well
                for (uint dd = 0; dd < D; dd++)
                {
                    string name = eqn->get_name();
                    switch(dd) {
                        case 0: name.append("_x"); break;
                        case 1: name.append("_y"); break;
                        case 2: name.append("_z"); break;
                        default: NEGF_EXCEPTION("strange dimension"); break;
                    }
                    vector<double> * component = new vector<double>;
                    component->resize(grid->get_num_vertices(), 0.0);
                    created_fields.push_back(component);
                    for (uint jj = 0; jj < grid->get_num_vertices(); jj++)
                        (*component)[jj] = eqn->get_value(D*jj+dd);
                    num_datasets++;
                    value_array.push_back(component);
                    datanames.push_back(name);
                    types.push_back(dat_unittypes[ii]);
                    locations.push_back("vertex");
                }
            } else {
                NEGF_EXCEPTION("Equation must be defined on vertices, edges or elements. Other possibility: Eqn is for a different grid.");
            }
        }
    }

    // write to .dat-file
    this->parser->write_dfise_dat(
        (this->resultfilename + ".dat").c_str(),
        this->grid,
        value_array,
        num_datasets,
        datanames,
        locations,
        types         );

    // clean up the created data vectors
    for (uint ii = 0; ii < created_fields.size(); ii++)
        delete created_fields[ii];
);}
