envName = "IAATEnv"

echo "****************************************************************************"
echo "*               Initial Allignemnt Assistant Tool - Installer              *"
echo "*               For assistance contact:                                    *"
echo "*                 - Andre Guimaraes: aguima1@lsu.edu                       *"
echo "*                 - _______________: _______________                       *"
echo "****************************************************************************"
echo ""
echo ""

echo "Welcome to the installer for the Initial Allignment Assistant Tool."

echo "Please provide the name of the preferred conda environment to be craeted (Press [Enter] for default: IAAT-Env)"

read envName

if [ "$envName" == "" ]; then
    envName="IAAT-Env"
fi

echo "Environment Name: $envName."

echo "Allowing conda on the shell."
source /cvmfs/ligo-containers.opensciencegrid.org/lscsoft/conda/latest/etc/profile.d/conda.sh

echo "Creating conda environment $envName (Might take a several minutes)."
conda env create --name $envName -f environmentIAAT.yml

echo "Activating the conda environment $envName."
conda activate $envName

echo "Installing ipykernel."
conda install -c anaconda ipykernel

echo "Adding $envName to the list of jupyter kernels."
python -m ipykernel install --user --name $envName --display-name "$envName"


