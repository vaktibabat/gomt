package main

import (
	"bufio"
	"fmt"
	"github.com/vaktibabat/gomt"
	"os"
	"strconv"
)

func check(err error) {
	if err != nil {
		panic(err)
	}
}

func main() {
	var my_state [624]uint32
	state_idx := 0
	leak_path := os.Args[1]

	fmt.Println("Guessing Game Solver")
	fmt.Println("Specify the path of the leak file as a command line arg")

	f, err := os.Open(leak_path)

	check(err)

	// Initialize a scanner to read the file line-by-line
	fScanner := bufio.NewScanner(f)
	fScanner.Split(bufio.ScanLines)

	my_state[0] = 0

	// Each line in the file is a leak of a number generated by the PRNG
	for fScanner.Scan() {
		// Convert to a number
		leak, err := strconv.Atoi(fScanner.Text())
		check(err)
		// Crack the corresponding state element
		my_state[state_idx] = mtlib.RestoreState(uint32(leak))

		state_idx += 1
	}

	my_mt := mtlib.MtFromState(my_state)

	for i := 0; i < 10; i++ {
		fmt.Printf("Next number: %d\n", my_mt.GenNext()%1000000)
	}
}
