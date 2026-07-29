/**
 * @module sage/groups
 * @description Group-theoretic algorithms and structures
 */

export {
  // Types
  type OperationType,
  type GroupElement,
  type MultiplicativeGroupElement,
  type AdditiveGroupElement,
  type GroupOps,
  // Constants
  MULTIPLICATION_NAMES,
  ADDITION_NAMES,
  // Functions
  parseGroupOps,
  multiple,
  bsgs,
  pohlig_hellman,
  discrete_log,
  discrete_log_generic,
  discrete_log_lambda,
  order_from_multiple,
  order_from_bounds,
  multiple_of_order,
  has_order,
  multiples,
  discrete_log_rho,
} from './generic.js';
