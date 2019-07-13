import asyncio


async def async_sh_job(cmd, sema):
    with (await sema):
        p = await asyncio.create_subprocess_shell(
            cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT)
        return (await p.communicate())[0].splitlines()


def async_batch_sh_jobs(cmd_list, thread=2):
    semaphore = asyncio.Semaphore(thread)
    if cmd_list:
        loop = asyncio.get_event_loop()
        coro_list = [async_sh_job(cmd, semaphore) for cmd in cmd_list]
        loop.run_until_complete(asyncio.wait(coro_list))
        loop.close()
