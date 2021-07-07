# -- coding:utf-8 --

import sys
import os
import glob
from subprocess import Popen, PIPE


def generate_args(binary, *params):
    arguments = [binary]
    arguments.extend(list(params))
    return arguments


def execute_binary(args):
    process = Popen(' '.join(args), shell=True, stdout=PIPE, stderr=PIPE)
    (std_output, std_error) = process.communicate()
    process.wait()
    rc = process.returncode

    return rc, std_output, std_error


def check_correctness(binary_path, data_graph_path, query_folder_path, expected_results):
    # find all query graphs.
    query_graph_path_list = glob.glob('{0}/*.graph'.format(query_folder_path))
    #print(query_graph_path_list)
    # check the correctness of Filtering.
    filter_type_list = ['DPiso','GQL']

    for filter_type in filter_type_list:
        for query_graph_path in query_graph_path_list:
            execution_args = generate_args(binary_path, '-d', data_graph_path, '-q', query_graph_path, '-filter', filter_type,
                                       '-order', 'DPiso', '-engine', 'DPiso')

            (rc, std_output, std_error) = execute_binary(execution_args)
            query_graph_name = os.path.splitext(os.path.basename(query_graph_path))[0]
            expected_embedding_num = expected_results[query_graph_name]

            if rc == 0:
                embedding_num = 0
                std_output = std_output.decode(encoding='utf-8',errors='strict')
                std_output_list = std_output.split('\n')
                #print(std_output_list)
                for line in std_output_list:
                    if '#Embeddings' in line:
                        embedding_num = int(line.split(':')[1].strip())
                        break
                if embedding_num != expected_embedding_num:
                    print('Filter type {0} {1} is wrong. Expected {2}, Output {3}'.format(filter_type, query_graph_name, expected_embedding_num,
                                                                      embedding_num))
                    exit(-1)
            else:
                print(rc)
                print('Filter type {0} {1} error.'.format(filter_type, query_graph_name))
                exit(-1)

        print('Filter type {0} passes the correctness check.'.format(filter_type))

    # check the correctness of ordering.
    order_type_list = ['DPiso','GQL']
    for order_type in order_type_list:
        for query_graph_path in query_graph_path_list:
            execution_args = generate_args(binary_path, '-d', data_graph_path, '-q', query_graph_path, '-filter', 'LDF',
                                       '-order', order_type, '-engine', 'LFTJ')

            (rc, std_output, std_error) = execute_binary(execution_args)
            query_graph_name = os.path.splitext(os.path.basename(query_graph_path))[0]
            expected_embedding_num = expected_results[query_graph_name]

            if rc == 0:
                embedding_num = 0
                std_output = std_output.decode(encoding='utf-8',errors='strict')
                std_output_list = std_output.split('\n')
                for line in std_output_list:
                    #找到报告embedding数目的输出行位置
                    if '#Embeddings' in line: 
                        embedding_num = int(line.split(':')[1].strip())
                        break
                if embedding_num != expected_embedding_num:
                    print('Order type {0} {1} is wrong. Expected {2}, Output {3}'.format(order_type, query_graph_name,
                            expected_embedding_num, embedding_num))
                    exit(-1)
            else:
                print('Order type {0} {1} error.'.format(order_type, query_graph_name))
                exit(-1)

        print('Order type {0} passes the correctness check.'.format(order_type))

    # check the correctness of ordering.
    order_type_list = ['DPiso','GQL']
    for order_type in order_type_list:
        for query_graph_path in query_graph_path_list:
            execution_args = generate_args(binary_path, '-d', data_graph_path, '-q', query_graph_path, '-filter', 'LDF',
                                       '-order', 'DPiso', '-engine', order_type)

            (rc, std_output, std_error) = execute_binary(execution_args)
            query_graph_name = os.path.splitext(os.path.basename(query_graph_path))[0]
            expected_embedding_num = expected_results[query_graph_name]

            if rc == 0:
                embedding_num = 0
                std_output = std_output.decode(encoding='utf-8',errors='strict')
                std_output_list = std_output.split('\n')
                for line in std_output_list:
                    #找到报告embedding数目的输出行位置
                    if '#Embeddings' in line: 
                        embedding_num = int(line.split(':')[1].strip())
                        break
                if embedding_num != expected_embedding_num:
                    print('Engine type {0} {1} is wrong. Expected {2}, Output {3}'.format(order_type, query_graph_name,
                            expected_embedding_num, embedding_num))
                    exit(-1)
            else:
                print('Engine type {0} {1} error.'.format(order_type, query_graph_name))
                exit(-1)

        print('Engine type {0} passes the correctness check.'.format(order_type))

if __name__ == '__main__':
    input_binary_path = './main'
    if not os.path.isfile(input_binary_path):
        print('The binary {0} does not exist.'.format(input_binary_path))
        exit(-1)

    # load expected results.
    input_expected_results = {}
    input_expected_results_file = '/workspace/match/test/expected_output.res'

    with open(input_expected_results_file, 'r') as f:
        for line in f:
            if line:
                result_item = line.split(':')
                input_expected_results[result_item[0].strip()] = int(result_item[1].strip())

    dir_path = os.path.dirname(os.path.realpath(__file__))
    input_data_graph_path = '{0}/data_graph/HPRD.graph'.format(dir_path)
    input_query_graph_folder_path = '{0}/query_graph/'.format(dir_path)

    input_data_graph_path = '/workspace/match/test/data_graph/HPRD.graph'
    input_query_graph_folder_path = '/workspace/match/test/query_graph/'
    
    check_correctness(input_binary_path, input_data_graph_path, input_query_graph_folder_path, input_expected_results)

