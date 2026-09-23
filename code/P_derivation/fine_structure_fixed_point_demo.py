#!/usr/bin/env python3
# Run this file to print the OPH fine-structure fixed-point value with a PDG/CERN hadron input.

from __future__ import annotations  # Keep annotations lightweight for this compact OPH demo.

import argparse  # Parse optional comparison inputs without putting target values in the solver.
from decimal import Decimal, localcontext  # Use decimal arithmetic so the displayed wave calculation is reproducible.


PRECISION = 90  # Use enough digits to expose the OPH fixed-point arithmetic.
ONE = Decimal(1)  # Represent the unit amplitude in reciprocal electromagnetic coupling.
TWO = Decimal(2)  # Represent the factor in the golden-ratio entropy equation.
FIVE = Decimal(5)  # Represent the discriminant of x^2 - x - 1 = 0.
PI = Decimal("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899")  # Provide pi for the Gaussian boundary wave width.
SOURCE_ANCHOR_INV_ALPHA_MZ = Decimal("128.30796547328624820996110874175671618724547618036535646005342169635117784168285644")  # Use the OPH source-chain electroweak anchor a0(P*); non-certified run output, independently recomputed to 27 significant digits.
LEPTON_TRANSPORT_DELTA_INV_ALPHA = Decimal("4.3093978664522040271317438975344894018487156605576773194711528089665680313257906466129")  # Use the OPH one-loop charged-lepton transport packet; non-certified run output, independently recomputed to 25 significant digits (convergence-budget limited, not certified).
PDG_DELTA_ALPHA_HAD_5_MZ = Decimal("0.02761")  # Declared input: the hadronic vacuum-polarization row Delta alpha_had^(5)(M_Z) = 0.02761 of the PDG electroweak review by Erler and Ferro-Hernandez, taken as printed; the PDG global fit carries its own hadronic input, so the row is varied within its stated uncertainty when the endpoint's sensitivity is checked.


def main() -> int:  # Keep the script runnable from the command line.
    parser = argparse.ArgumentParser(description="Print the OPH alpha fixed point using a PDG hadron input.")  # Define the demo command line.
    parser.add_argument("--compare-alpha-inv", default=None, help="Optional inverse-alpha comparison value; never enters the solve.")  # Keep target values outside the computation path.
    args = parser.parse_args()  # Read optional comparison arguments from the user.
    with localcontext() as ctx:  # Isolate the high-precision OPH fixed-point arithmetic.
        ctx.prec = PRECISION  # Set precision for phi, sqrt(pi), transport, P, and alpha.
        phi = +((ONE + FIVE.sqrt()) / TWO)  # Compute the OPH golden-ratio screen balance.
        sqrt_pi = +(PI.sqrt())  # Compute the Gaussian wave-width normalization on the boundary.
        source_plus_leptons = +(SOURCE_ANCHOR_INV_ALPHA_MZ + LEPTON_TRANSPORT_DELTA_INV_ALPHA)  # Combine OPH electroweak anchor with OPH leptonic transport.
        alpha_inv_endpoint = +(source_plus_leptons / (ONE - PDG_DELTA_ALPHA_HAD_5_MZ))  # Convert the PDG hadronic vacuum polarization into the Thomson endpoint convention.
        alpha_from_endpoint = +(ONE / alpha_inv_endpoint)  # Convert the Thomson endpoint width into alpha.
        pixel = +(phi + sqrt_pi / alpha_inv_endpoint)  # Solve the OPH fixed point P = phi + sqrt(pi) / A_Th(P).
        alpha_from_pixel = +((pixel - phi) / sqrt_pi)  # Read the same coupling from the pixel detuning.
        alpha_inv_from_pixel = +(ONE / alpha_from_pixel)  # Convert the geometric alpha back to inverse alpha.
        fixed_point_residual = +(pixel - (phi + sqrt_pi / alpha_inv_endpoint))  # Check the OPH pixel fixed-point equation.
        coupling_residual = +(alpha_from_pixel - alpha_from_endpoint)  # Check that inside wave width equals outside detuning.
        print("OPH fine-structure fixed-point demo with PDG/CERN hadron input")  # Label the run as the no-target-input calculation.
        print(f"phi                                  = {phi}")  # Print the golden-ratio OPH screen balance.
        print(f"sqrt(pi)                             = {sqrt_pi}")  # Print the Gaussian boundary wave factor.
        print(f"OPH source anchor a0(P*)              = {SOURCE_ANCHOR_INV_ALPHA_MZ}")  # Print the OPH electroweak source anchor.
        print(f"OPH lepton transport delta            = {LEPTON_TRANSPORT_DELTA_INV_ALPHA}")  # Print the source lepton transport packet.
        print(f"PDG Delta alpha_had^(5)(M_Z)          = {PDG_DELTA_ALPHA_HAD_5_MZ}")  # Print the empirical hadron input used by the solve.
        print(f"A_Th(P*) = alpha^-1(0)                = {alpha_inv_endpoint}")  # Print the resulting inverse fine-structure constant.
        print(f"alpha(0)                             = {alpha_from_pixel}")  # Print the fine-structure constant from the pixel.
        print(f"P*                                   = {pixel}")  # Print the pixel area fixed point.
        print(f"alpha^-1 from P check                = {alpha_inv_from_pixel}")  # Print the reciprocal recovered from the pixel.
        print(f"P fixed-point residual               = {fixed_point_residual}")  # Print closure error for P.
        print(f"alpha matching residual              = {coupling_residual}")  # Print closure error for alpha.
        if args.compare_alpha_inv is not None:  # Run comparison math only when the user supplies a target value.
            compare_alpha_inv = +Decimal(args.compare_alpha_inv)  # Parse the comparison target after the OPH solve is finished.
            compare_gap = +(alpha_inv_endpoint - compare_alpha_inv)  # Compare the no-target output against the supplied target.
            required_delta_had = +(ONE - source_plus_leptons / compare_alpha_inv)  # Compute the hadron delta needed in this convention.
            required_pre_hadronic = +(compare_alpha_inv * (ONE - PDG_DELTA_ALPHA_HAD_5_MZ))  # Compute the pre-hadronic inverse-alpha value implied by the target and PDG hadrons.
            required_source_anchor = +(required_pre_hadronic - LEPTON_TRANSPORT_DELTA_INV_ALPHA)  # Remove the OPH lepton packet to get the source anchor implied by the target.
            anchor_gap = +(required_source_anchor - SOURCE_ANCHOR_INV_ALPHA_MZ)  # Compute the missing same-scheme electroweak anchor bridge.
            print(f"comparison alpha^-1 target           = {compare_alpha_inv}")  # Print the external comparison target.
            print(f"comparison gap                       = {compare_gap}")  # Print output minus comparison target.
            print(f"required Delta alpha_had in this convention = {required_delta_had}")  # Print the hadronic delta implied by the target.
            print(f"required pre-hadronic value with PDG hadrons = {required_pre_hadronic}")  # Print the source-plus-lepton value implied by the target.
            print(f"required source anchor with PDG hadrons     = {required_source_anchor}")  # Print the anchor implied by the target and PDG hadrons.
            print(f"same-scheme source-anchor gap              = {anchor_gap}")  # Print the electroweak scheme bridge gap.
    return 0  # Signal that the OPH fixed-point calculation completed.


if __name__ == "__main__":  # Execute only when run as a script.
    raise SystemExit(main())  # Return the command-line exit code.
