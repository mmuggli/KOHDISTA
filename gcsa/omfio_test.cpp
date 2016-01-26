#include "omfio.h"

int main(int argc, char** argv)
{
    Omfio rmaps(argv[1]);
    rmaps.dump();

    return 0;
}
