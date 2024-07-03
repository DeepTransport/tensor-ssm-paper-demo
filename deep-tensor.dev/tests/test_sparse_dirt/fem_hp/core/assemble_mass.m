function mass = assemble_mass(mesh, local_elem, fill_i, f, fill_d)
%assemble_mass
%
% assemble the mass matrix 
%
% Tiangang Cui, 17/Oct/2016

loc_f = local_elem.f_quad_pts*f(mesh.node_map');
fill  = (local_elem.mass_at_q.*local_elem.quad_weights(:)')*(loc_f.*mesh.detJ(:)');
mass  = accumarray(fill_i, fill(:), fill_d, [], [], true);

end

%function mass = assemble_mass(FEM, f_at_q)
