#create unitttest for ims3d.py
import unittest
import os
import tempfile
import sys

#import ims3d and ims3d_harv from .. directory
#  we assume we are in IMS3DPATH/tests directory
IMS3DPATH = os.path.dirname(os.path.abspath(__file__)) + '/../'
sys.path.append(IMS3DPATH)
import ims3d
import ims3d_harv

class TestIMS3D(unittest.TestCase):
    def test_IMS3D_exists(self):
        result = ims3d.main
        self.assertIsNotNone(result)

    def test_IMS3D_tutorial_com_file(self):
        reference_com_file = os.path.join(IMS3DPATH, 'tests', 'references', 'input_batch_00000_ref.com')
        generated_com_file = 'input_batch_00000.com'
        # append the current directory to a filename
        sys.argv = ['ims3d.py', '-r', '1', os.path.join(IMS3DPATH, 'tutorial', 'naphtalene.xyz')]
        #create a temporary directory and move into it
        with tempfile.TemporaryDirectory() as tempdir:
            os.chdir(tempdir)
            # Call the main function of ims3d
            result = ims3d.main()
            # Check if the generated file is similar to the reference file
            generated_content, reference_content = generate_and_compare_contents(reference_com_file, generated_com_file)
            self.assertEqual(generated_content, reference_content)

def generate_and_compare_contents(reference_com_file, generated_com_file):
    with open(generated_com_file, 'r') as generated_file:
        generated_content = generated_file.read()
    with open(reference_com_file, 'r') as reference_file:
        reference_content = reference_file.read()
        # if contents differ, print the differences
    if generated_content != reference_content:
        print("Generated content does not match reference content.")
        print("Generated content:")
        print(generated_content)
        print("Reference content:")
        print(reference_content)
    return generated_content,reference_content

class TestIMS3DHarv(unittest.TestCase):
    def test_IMS3D_harv_exists(self):
        # Example test case, replace with actual test logic
        result = ims3d_harv.main
        self.assertIsNotNone(result)  # Check if the result is not None
    
    def test_IMS3D_tutorial_log_file(self):
        reference_log_file = os.path.join(IMS3DPATH, 'tests', 'references', 'input_batch_00000_ref.log')
        generated_dat_file = 'ims.dat'
        reference_dat_file = os.path.join(IMS3DPATH, 'tests', 'references', 'ims_ref.dat')
        sys.argv = ['ims3d_harv.py', reference_log_file]
        #create a temporary directory and move into it
        with tempfile.TemporaryDirectory() as tempdir:
            os.chdir(tempdir)
            # Call the main function of ims3d
            result = ims3d_harv.main()
            # Check if the generated file is similar to the reference file
            generated_content, reference_content = generate_and_compare_contents(reference_dat_file, generated_dat_file)
            self.assertEqual(generated_content, reference_content)

if __name__ == '__main__':
    unittest.main()
