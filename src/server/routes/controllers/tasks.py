# stdlib
import logging

# third party lib
from flask import request

# internal lib
from lib import tasks
from lib import networking
from lib.settings import Commands, Options, Thresholds
from server.lib import task_init, task_qc, task_pca, task_ass
from lib.client_registry import Registry

# TODO helper for parameter validation and default specification


def list_tasks():
    tsks = tasks.task_list
    return networking.create_response(200, tsks)


def start_task(task_name):
    logging.info(f'Got command to start {task_name}, starting...')
    args = {}

    for key, val in request.json.items():
        args[key.upper()] = val

    if task_name == Commands.INIT:
        task_init.start_init_task()

    elif task_name.startswith(Commands.QC):
        args[Options.HWE] = args.get("QC_HWE", Thresholds.QC_hwe)
        args[Options.MAF] = args.get("QC_MAF", Thresholds.QC_maf)
        logging.info(f"Specified Filters :{args}")
        task_qc.start_client_qc_task(args)
        task_qc.start_local_qc_task(args)

    elif task_name.startswith(Commands.PCA):
        args["PCA_PCS"] = int(args.get("PCA_PCS", Thresholds.PCA_pcs))
        task_pca.Position_reporter.get_instance(args)
        if not task_pca.ready_to_decompose():
            if not task_pca.filtered():
                args[Options.MAF] = args.get("PCA_MAF", Thresholds.PCA_maf)
                if "PCA_LD_WINDOW" not in args:  # default parameters
                    args["PCA_LD_WINDOW"] = Thresholds.PCA_ld_window
                if "PCA_LD_THRESHOLD" not in args:  # default parameters
                    args["PCA_LD_THRESHOLD"] = Thresholds.PCA_ld_threshold
                args[Options.LD] = [args["PCA_LD_WINDOW"], args["PCA_LD_THRESHOLD"]]
                logging.info(f"Specified pruning filters :{args}")
                task_pca.start_pca_filters(args)
            else:
                logging.info("Reporting Filtered Sites")
                task_pca.Position_reporter.get_instance().report_pos()
        else:
            logging.info("starting eigen decomposition")
            task_pca.eigenDecompose(n_components=args["PCA_PCS"])

    elif task_name == Commands.ASSO:
        args["ASSO_PCS"] = args.get("ASSO_PCS", Thresholds.ASSO_pcs)
        logging.info("Starting Associations")
        task_ass.LogisticAdmm.get_instance(args, active=2)
    
    elif task_name == Commands.ECHO:
        counts = args.get("ECHO_COUNTS", Thresholds.ECHO_COUNTS)
        echo = task_init.Echo.get_instance(counts)

    return networking.create_response(200, f'Started task {task_name}')


def start_subtask(task_name, subtask_name, client_name):
    if task_name == Commands.INIT:
        if subtask_name == 'POS':
            task_init.store_positions(request.data, client_name)
        elif subtask_name == 'COUNT':
            task_init.store_counts(request.data, client_name)

    elif task_name.startswith(Commands.QC):
        if subtask_name == "FIN":
            if task_qc.filter_finished(client_name, Commands.QC):
                logging.info("Done with QC.")

    elif task_name.startswith(Commands.PCA):
        if subtask_name == "FIN":
            if task_qc.filter_finished(client_name, Commands.PCA):
                logging.info("Done with PCA filters. Initiating pruning")
                reset_states("PRUNE")
                ld_agg = task_pca.CovarianceAggregator.get_instance(len(Registry.get_instance().list_clients()), 50)
                # send message to start LD pruning
                ld_agg.send_request({})
        elif subtask_name == "LD":
            ld_agg = task_pca.CovarianceAggregator.get_instance(len(Registry.get_instance().list_clients()), 50)
            ld_agg.update(request.data)
        elif subtask_name == "PCAPOS":
            task_pca.Position_reporter.get_instance().report_pos()
        elif subtask_name == "COV":
            task_pca.store_covariance(client_name, request.data)

    elif task_name.startswith(Commands.ASSO):
        ass_agg = task_ass.LogisticAdmm.get_instance({}, active=2)
        if subtask_name == "adjust":
            ass_agg.update_stats(request.data)
        elif subtask_name == "estimate":
            ass_agg.update(request.data)
        elif subtask_name == "pval":
            ass_agg.update_pval(request.data)
        elif subtask_name == "hessians":
            model, have_all_info = ass_agg.newton_stats_update(request.data)
            if have_all_info:
                ass_agg.newton_iter(model)
        elif subtask_name == "valback":
            ass_agg.collect_likelihoods(request.data)
    elif task_name == Commands.ECHO:
        if subtask_name == "ITR":
            echo = task_init.Echo.get_instance(1)
            echo.echo(client_name)
            avg_t = echo.echo(client_name)
            if avg_t is not None:
                logging.info(f"Avg echo time={avg_t}")

    return networking.create_response(200)


def reset_states(state):
    instance = Registry.get_instance()
    for client in instance.list_clients():
        instance.set_client_state(client["name"], state)


def next_task():
    tsk = tasks.TaskReg.get_instance()
    msg = tsk.get_up_task()
    return networking.create_response(200, msg)
