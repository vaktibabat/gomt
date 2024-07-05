package mtlib

// Word size (i.e. the size of the numbers we generate in bits)
const w uint = 32

// Degree of reccurence
const n int = 624

// Offset used in the reccurence
const m uint = 397

// The last row of the matrix A
const a uint = 0x9908B0DF

// Constant for initialization
const f uint = 1812433253

// Bitmask for taking the upper w-r=1 bits
const upper_mask uint = 0x80000000

// Bitmask for taking the lower r bits
const lower_mask uint = 0x7FFFFFFF

// Constants for the tempering matrix T
const u uint = 11
const s uint = 7
const b uint = 0x9D2C5680
const t uint = 15
const c uint = 0xEFC60000
const l uint = 18

// The current state of the Mersenne Twister
type mt_state struct {
	state     [n]uint32
	state_idx uint
}

// Initialize a new MT
func NewMt(seed uint32) *mt_state {
	var state [n]uint32
	mt := mt_state{state: state, state_idx: 0}
	mt.init_state(seed)

	return &mt
}

// Initialize a new MT given some state
func MtFromState(state [n]uint32) *mt_state {
	mt := mt_state{state: state, state_idx: 0}

	return &mt
}

// Initialize the state of the MT according to some seed
func (mt *mt_state) init_state(seed uint32) {
	//  First element of the initial state is the seed
	mt.state[0] = seed

	// All of the next elements are initialized
	// with a reccurence
	for i := 1; i < n; i++ {
		prev := uint(mt.state[i-1])

		mt.state[i] = uint32(f*(prev^(prev>>(w-2)))) + uint32(i)
	}

	mt.state_idx = 0
}

// Generate a random number from the MT and mutate the oldest number in the state
func (mt *mt_state) GenNext() uint32 {
	// The current state index
	i := int(mt.state_idx)
	// Compute the concatenation between the upper 1 bit of state[i] and lower 31 bits of state[i + 1]
	y := (mt.state[i] & uint32(upper_mask)) | (mt.state[(i+1)%n] & uint32(lower_mask))
	// Mutate the state by adding the i+m % n th number in the state to the matrix product
	y_lsb := y & 1
	var mat_prod uint

	// If the LSB of y is 0, the last row of A, which is the bits of a, won't be added to the result
	if y_lsb == 0 {
		mat_prod = 0
	} else {
		// Otherwise, we multiply a 1 with every bit of a, which yields a itself
		mat_prod = a
	}
	// If the LSB of y is 0, it will not be added by
	mt.state[i] = mt.state[(i+int(m))%n] ^ (y >> 1) ^ uint32(mat_prod)
	// Get the output by computing the product of the state with the tempering matrix
	out := mt.state[i]
	out ^= out >> uint32(u)
	out ^= (out << s) & uint32(b)
	out ^= (out << t) & uint32(c)
	out ^= out >> uint32(l)
	// Update the state index
	mt.state_idx = uint((i + 1) % n)

	return out
}

// Given some output from the PRNG, restore the corresponding element in the state
// e.g. given output 5, we can restore element 5 (or 4 if we're using zero-based indexing)
func RestoreState(out uint32) uint32 {
	// MT generates the output by applying an *invertible* tempering
	// transformation to the state element
	tempered_state := out
	// Inverse of "out ^= out >> uint32(l)"
	tempered_state ^= tempered_state >> uint32(l)
	// Inverse of "out ^= (out << t) & uint32(c)"
	tempered_state ^= (tempered_state << t) & uint32(c)
	// Inverse of "out ^= (out << s) & uint32(b)"
	tempered_state ^= ((tempered_state << 28) & 0x10000000) ^ ((tempered_state << 21) & 0x14200000) ^ ((tempered_state << 14) & 0x94284000) ^ ((tempered_state << s) & uint32(b))
	// Inverse of "out ^= out >> uint32(u)"
	original_state := tempered_state ^ (tempered_state >> 11) ^ (tempered_state >> 22)

	return original_state
}
