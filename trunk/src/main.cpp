#include "smoke3D.h"
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char * argv[]) {
	smoke3D::init();
	while(1) smoke3D::simulateStep();
	return 0;
}
