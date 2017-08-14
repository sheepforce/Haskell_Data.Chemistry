-- test the energy_calculator with the Müller-Brown Potential
let educt = BLAS.vector [0.6211317901475866, 2.823674628779631e-2]
let prod = BLAS.vector [-0.5591208056179865, 1.4407836667987057]

let molecule_educt = QcInMolecule {headeropts = "", name = "educt", memory = (15, "Gb"), atoms = ["A", "B"], geom_cart = educt, multiplicity = 1, charge = 0, inputopts = "energy(\"scf\")"}
let molecule_prod = QcInMolecule {headeropts = "", name = "product", memory = (15, "Gb"), atoms = ["A", "B"], geom_cart = prod, multiplicity = 1, charge = 0, inputopts = "energy(\"scf\")"}
                                  
energy_calculator 0 molecule_educt Test educt
gradient_calculator 0 molecule_educt Test educt

-- test the Psi4 energy and gradient wrapper calculator
input <- T.readFile "/home/phillip/Ethen.xyz"
let ethen_xyz = fromRight $ parseOnly XYZ.xyzParser input
let ethen_vec = XYZ.coord2Vec ethen_xyz
let ethen_qc = QcInMolecule {headeropts = "", name = "Ethen", memory = (15, "Gb"), atoms = XYZ.getElements ethen_xyz , geom_cart = ethen_vec , multiplicity = 1, charge = 0, inputopts = "set {\n basis def2-svp\n}\ngradient(\"scf\")"}

energy_calculator ethen_qc Psi4 "EthenGrad" "/usr/bin/psi4" 4 ethen_vec
gradient_calculator ethen_qc Psi4 "EthenGrad" "/usr/bin/psi4" 4 ethen_vec 

-- test the interpolation of images in the ridge algorithm with SN substitution
edukt_text <- T.readFile "/home/phillip/Edukt.xyz"
produkt_text <- T.readFile "/home/phillip/Produkt.xyz"
let educt_xyz = fromRight $ parseOnly XYZ.xyzParser edukt_text
let product_xyz = fromRight $ parseOnly XYZ.xyzParser produkt_text
let elements = XYZ.getElements educt_xyz
let educt_vec = XYZ.coord2Vec educt_xyz
let product_vec = XYZ.coord2Vec product_xyz
template_text <- T.readFile "/home/phillip/Template.psi"
let template = T.unpack template_text

--ridge_optimization elements educt_vec product_vec 1 (-1) 5 2 50 0.05 0.05 0.1 1.0e-6 0.02 (1.0e-4, 5.0e-4) (Psi4, template) BFGS2 4 High 0 [] []
ridge_optimization elements educt_vec product_vec 1 (-1) 5 10 50 0.25 1.0 0.1 1.0e-6 0.02 (1.0e-4, 5.0e-4) (Psi4, template) BFGS2 4 High 0 [] []


-- test the matrix printing function
lAK = BLAS.vector [-200, -100, -170, 15]
ak = BLAS.vector [-1, -1, -6.5, 0.7]
bk = BLAS.vector [0, 0, 11, 0.6]
ck = BLAS.vector [-10, -10, -6.5, 0.7]
x0k = BLAS.vector [1, 0, -0.5, -1]
y0k = BLAS.vector [0, 0.5, 1.5, 1]

let start = BLAS.vector [0.6211317901475866,  2.823674628779631e-2]
let end = BLAS.vector [-0.5591208056179865, 1.4407836667987057]
let reaction_vec = end - start
let solMat = steepestDescent 1.0e-6 1000 0.02 1.0e-5 (negate . muellerBrown (lAK, ak, bk, ck, x0k, y0k)) (negate . (vectorProjection reaction_vec) . muellerBrown' (lAK, ak, bk, ck, x0k, y0k)) end
