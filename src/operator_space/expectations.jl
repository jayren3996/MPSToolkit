"""
    pauli_trace(rho)

Return the physical trace of a vectorized operator stored as a Pauli-basis `MPS`.

In the normalized Pauli convention `tr(ρ) = (√2)^N c_{I…I}`, where `c_{I…I}` is the
amplitude of the all-identity string, so the trace is one identity-environment contraction.

# Arguments
- `rho`: Operator-space `MPS` on dimension-4 Pauli sites.

# Returns
- The trace as a `ComplexF64` (real up to numerical noise for Hermitian operators).
"""
function pauli_trace(rho::MPS)
  _validate_pauli_operator_space(rho, 1, length(rho))
  return 2.0^(length(rho) / 2) * scalar(_right_identity_environment(rho, 1))
end
